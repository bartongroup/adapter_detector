import os
import numpy as np

import pysam
import click

from .utils import get_fast5_read_id_mapping, get_output_bam_filenames
from .model import load_adapter_model
from .signal import get_fast5_fiveprime


def process_batch(batch, reads, 
                  model, threshold,
                  batch_size, squiggle_size):
    '''
    Make a prediction for a batch of signals
    '''
    batch = np.array(batch).reshape(batch_size, squiggle_size, 1)
    preds = model.predict(batch).squeeze()
    for score, read in zip(preds, reads):
        read.set_tag('ac', '{:.3f}'.format(score), value_type='f')
        yield score > threshold, read


def bam_filter(bam_fn, read_id_fast5_filemap, model,
               pass_output_bam_fn, fail_output_bam_fn, threshold=0.5,
               batch_size=5000, squiggle_size=2000):
    '''
    Filter a bam file using a model trained to detect adapters in
    signal from 5' end.
    
    Parameters:
    ----------
    
        bam_fn: str, required
          Path to input bam file

        read_id_fast5_filemap: dict, required
          Dict of read_id: fast5 filepath key value pairs

        model: keras.models.Model, required
          Model used to make predictions. Must accept input of
          shape (batch_size, squiggle_size, 1)

        pass_output_bam_fn: str, required
          Bam file path to write reads passing filter to.
          Will be overwritten

        fail_output_bam_fn: str, required
          Bam file path to write reads failing filter to.
          Will be overwritten

        threshold: float, optional, default: 0.5
          Threshold at which to cutoff pass/fail

        batch_size, int, optional, default: 5000
          Number of signals to pass to model prediction
          at a time

        squiggle_size: int, optional, default: 2000
          Length of signal to take from 5' end of RNA fast5 signal,
          to pass to model.

    Returns:
    --------
    
        None
    '''
    in_bam = pysam.AlignmentFile(bam_fn, mode='rb')
    # create two output bams
    out_bams = {
        True: pysam.AlignmentFile(pass_output_bam_fn,
                                  mode='wb', template=in_bam),
        False: pysam.AlignmentFile(fail_output_bam_fn,
                                   mode='wb', template=in_bam)
    }
    batch = []
    batch_reads = []
    for record in in_bam.fetch():
        read_id = record.query_name
        try:
            f5_fn = read_id_fast5_filemap[read_id]
        except KeyError:
            continue
        batch_reads.append(record)
        batch.append(get_fast5_fiveprime(f5_fn, squiggle_size))
        if len(batch) == batch_size:
            for passes, read in process_batch(batch, batch_reads,
                                              model, threshold,
                                              batch_size, squiggle_size):
                out_bams[passes].write(read)
            batch = []
            batch_reads = []
    if len(batch):
        for passes, read in process_batch(batch, batch_reads,
                                          model, threshold,
                                          batch_size, squiggle_size):
            out_bams[passes].write(read)
    in_bam.close()
    pos_bam.close()
    neg_bam.close()


@click.group()
def cli():
    pass

@cli.command()
@click.option('-b', '--bam-file', required=True, type=str)
@click.option('-s', '--summary-files-basedir', required=True, type=str)
@click.option('-f', '--fast5-files-basedir', required=True, type=str)
@click.option('-o', '--output-bam-basename', required=True, type=str)
@click.option('-m', '--model-file', required=False, default=None)
@click.option('-t', '--threshold', required=False, default=0.5, type=float)
@click.option('-n', '--batch-size', required=False, default=5000, type=int)
@click.option('-l', '--squiggle-size', required=False, default=2000, type=int)
def filter_bam(bam_file,
               summary_files_basedir,
               fast5_files_basedir,
               output_bam_basename,
               model_file,
               threshold,
               batch_size,
               squiggle_size):
    read_id_f5_filemap = get_fast5_read_id_mapping(
        summary_files_basedir,
        fast5_files_basedir,
    )
    pass_output_bam_fn, fail_output_bam_fn = get_output_bam_filenames(
        output_bam_basename
    )
    model = load_adapter_model(model_fn)
    bam_filter(
        bam_file,
        read_id_f5_filemap,
        model,
        pass_output_bam_fn,
        fail_output_bam_fn,
        threshold,
        batch_size,
        squiggle_size
    )
    

if __name__ == '__main__':
    run_bam_filter()