#! /usr/bin/env python

import os
import sys
import glob

from optparse import OptionParser

from myutils.BetterConfigParser import BetterConfigParser
from myutils.FileList import FileList
from myutils.FileLocator import FileLocator
from myutils.copytreePSI import filelist
from myutils.sample_parser import ParseInfo


def get_config(opts):
    # From submit.py

    pathconfig = BetterConfigParser()
    pathconfig.read('%sconfig/paths.ini' % opts.tag)

    try:
        _configs = [x for x in pathconfig.get('Configuration', 'List').split(" ") if len(x.strip()) > 0]
        if 'volatile.ini' in _configs:
            _configs.remove('volatile.ini')
        configs = ['%sconfig/' % (opts.tag) + c for c in _configs]

    except Exception as e:
        print("\x1b[31mERROR:" + str(e) + "\x1b[0m")
        print("\x1b[31mERROR: configuration file not found. Check config-tag specified with -T and presence of '[Configuration] List' in .ini files.\x1b[0m")
        raise Exception("ConfigNotFound")

    # read config
    config = BetterConfigParser()
    config.read(configs)

    return config


def add_to_config(condor_config, template, sample, config, locator):

    if os.environ.get('XBBDEBUG'):
        print('Adding %s:' % sample)

    sampledir = os.path.join(config.get('Directories', 'CONDORout'), sample)
    if not os.path.exists(sampledir):
        os.makedirs(sampledir)

    for part, infile in enumerate(filelist(config.get('Directories', 'samplefiles'), sample)):

        job = {
            'log': '%s_part%s' % (sample, part),
            'part': part,
            'sample': sample,
            'filelist': FileList.compress(infile),
            'outfile': locator.getFilenameAfterPrep(infile)
            }

        if os.path.exists(os.path.join(sampledir, job['outfile'])):
            continue

        condor_config.write(template.format(**job))


if __name__ == '__main__':

    parser = OptionParser()

    parser.add_option('-T', '--tag', dest='tag', default='default',
                      help='Tag to run the analysis with, example \'8TeV\' uses 8TeVconfig to run the analysis')
    parser.add_option('-S','--samples',dest='samples',default='*', help='samples you want to run on')
    parser.add_option('-o', '--output', dest='output', default='condor', help='output prefix')

    (opts, args) = parser.parse_args(sys.argv)

    config = get_config(opts)
    filelocator = FileLocator(config)
    parseinfo = ParseInfo(samples_path=config.get('Directories', 'PREPin'), config=config)

    with open('batch/condor/mit_header.sub', 'r') as header_file:
        header = header_file.read()

    logdir = os.path.join('/home/dabercro/public_html/xbb', config.get('Directories', 'Dname'))

    if not os.path.exists(logdir):
        os.makedirs(logdir)

    with open('batch/condor/template.sub', 'r') as template_file:
        template = template_file.read().format(
            logdir=logdir,
            tag=opts.tag,
            outdir=config.get('Directories', 'SYSout'),
            condorout=config.get('Directories', 'CONDORout'),
            log='{log}', part='{part}', sample='{sample}',
            filelist='{filelist}', outfile='{outfile}'
            )

    with open('%s_%s.cfg' % (opts.output, opts.tag), 'w') as condor_config:

        condor_config.write(header)

        for sample_file in glob.iglob('%s/%s.txt' % (config.get('Directories', 'samplefiles'), opts.samples)):

            if sample_file.endswith('.root.txt'):
                continue

            sample = os.path.basename(sample_file).split('.')[0]

            samples = parseinfo.find(sample)

            if os.environ.get('XBBDEBUG'):
                print(samples)

            if len(samples) == 1:
                add_to_config(condor_config, template, sample, config, filelocator)
