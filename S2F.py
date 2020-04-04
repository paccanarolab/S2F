"""

S2F - Main Script

This is the main entry point for S2F, all the commands can be run
by running this script using python 3.X

"""

import argparse
import commands
import json
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Sequence to Function'
    )
    subparsers = parser.add_subparsers(help='sub-command help',
                                       dest='subcommand')

    predict = subparsers.add_parser(
        'predict',
        description='S2F predict: given a fasta file, predicts ' +
                    'protein function using S2F',
        help='predict command'
    )
    predict.set_defaults(func=commands.predict)
    predict.add_argument('--run-config',
                         help='path to the run configuration file ' +
                              '(overrides all other arguments)',
                         default='arguments')
    predict.add_argument('--config-file',
                         help='location of the installation configuration ' +
                              'file that will be loaded. If not provided, '
                              'the default configuration file will be loaded',
                         default=os.path.join(os.path.dirname(
                                                  os.path.abspath(__file__)),
                                              's2f.conf'))
    predict.add_argument('--alias', help='Name of the prediction run')
    predict.add_argument('--obo', help='Path to the go.obo file to use')
    predict.add_argument('--fasta', help='Path to the protein sequencefile')
    predict.add_argument('--cpu',
                         help='Number of CPUs to use for ' +
                              'parallelisable computations', default='infer')
    predict.add_argument('--interpro-output',
                         help='manually provide InterPro output file and' +
                              'therefore avoid its computation',
                         default='compute')
    predict.add_argument('--hmmer-output',
                         help='manually provide HMMer output file and' +
                              ' therefore avoid its computation',
                         default='compute')
    predict.add_argument('--transfer-blacklist',
                         help='path to a file that containes a list of ' +
                              'identifiers from which no' +
                              ' links will be transferred. The file should ' +
                              'have one identifier per line',
                         default='')
    predict.add_argument('--hmmer-blacklist',
                         help='path to a file that contains a list of '
                              'identifiers that will be ignored ' +
                              'from the HMMer result list. The file should ' +
                              'have one identifier per line',
                         default='compute')
    predict.add_argument('--graph-collection',
                         help='provide a STRING graph collection manually ' +
                              'to avoid building one',
                         default='compute')
    predict.add_argument('--combined-graph',
                         help='manually provide a combined graph to avoid ' +
                              'its construction ' +
                              '(overwrites --graph-collection)',
                         default='compute')
    predict.add_argument('--homology-graph',
                         help='manually provide a homology graph file, ' +
                              'avoiding its computation',
                         default='compute')
    predict.add_argument('--goa-clamp',
                         help='provide a set of ground truth ' +
                              'functional associations, these ' +
                              'associations will be clamped (set ' +
                              'to 1) into the diffusion seed. The ' +
                              ' format of this file is: ' +
                              'PROTEIN_ID<tab>GO_ID',
                         default='compute')
    predict.add_argument('-ua', '--unattended',
                         help='The prediction will not be interactive, and '
                              'all configurations will be accepted without '
                              'prompting the user.',
                         action='store_true')

    # TODO: perhaps it is a good idea to add the BLAST executable
    # as a parameter as well
    install = subparsers.add_parser(
        'install',
        description='Installs the necessary requirements',
        help='install S2F and generate a configuration file'
    )
    install.set_defaults(func=commands.install)
    install.add_argument('--installation-directory',
                         help='Path to the installation directory for S2F.',
                         default='~/.S2F')
    install.add_argument('--config-file',
                         help='location of the configuration file ' +
                              'that will be created. If not provided,' +
                              'the default configuration file will be loaded',
                         default=os.path.join(
                              os.path.dirname(os.path.abspath(__file__)),
                              's2f.conf'))
    install.add_argument('--interpro',
                         help='manually provide the path to the ' +
                              '"iprscan" executable to avoid passing this ' +
                              'parameter to the other commands every time.',
                         default='iprscan')
    install.add_argument('--hmmer',
                         help='manually provide the route to the "phmmer"' +
                              'executable to avoid passing this '
                              'parameter to the other commands every time.',
                         default='phmmer')
    install.add_argument('--blastp',
                         help='manually provide the path to the ' +
                              '"blastp" command in the system. ' +
                              'If not provided, S2F will assume ' +
                              'that the executable is available system-wide.',
                         default='blastp')
    install.add_argument('--makeblastdb',
                         help='manually provide the path to the ' +
                              '"makeblastdb" command in the system. ' +
                              ' If not provided, S2F will assume that ' +
                              'the executable is available system-wide.',
                         default='makeblastdb')
    install.add_argument('--string-links',
                         help='manually provide the path to the STRING' +
                              'interactions DB, it must be the full path' +
                              'to either "protein.links.full.vX.x.txt.gz" ' +
                              'or "protein.links.detailed.vX.x.txt.gz". If '
                              'not provided, the installation script ' +
                              'will attempt to download the full '
                              'database using the wget command.',
                         default='download')
    install.add_argument('--string-sequences',
                         help='manually provide the path to the ' +
                              'STRING sequences database, it must ' +
                              ' be the full path to the ' +
                              ' "protein.sequences.vX.x.fa.gz" file. ' +
                              'If not provided, the installation script will' +
                              ' attempt to download it using "wget".',
                         default='download')
    install.add_argument('--string-species',
                         help='manually provide the path to the STRING' +
                              'species list, it must be the full path to the' +
                              ' "species.vX.x.txt" file. If not provided, ' +
                              ' the installation script will attempt to '
                              'download it using the wget command.',
                         default='download')
    install.add_argument('--uniprot-swissprot',
                         help='manually provide the path to the UniProt' +
                              'SwissProt sequences, it must be the full path' +
                              ' to the "goa_uniprot_all.gaf.gz" file. If not' +
                              'provided, the installation script will ' +
                              'attempt to download it using the wget command.',
                         default='download')
    install.add_argument('--uniprot-goa',
                         help='manually provide the path to the UniProt ' +
                              'GOA, it must be the full path to the ' +
                              '"goa_uniprot_all.gaf.gz" file. If not ' +
                              'provided, the installation script will ' +
                              'attempt to download it using the wget command.',
                         default='download')
    install.add_argument('--evidence-codes',
                         help='manually provide a list of evidence ' +
                              'codes that will be used to filter the UniProt' +
                              ' GOA. If not provided, S2F will be installed ' +
                              'using only experimental evidence codes.',
                         default='experimental')

    combine = subparsers.add_parser(
        'combine',
        description='S2F combine: run the network combination ' +
                    'on manually provided files',
        help='combine command'
    )
    combine.set_defaults(func=commands.combine)
    combine.add_argument('--fasta',
                         help='Path to the organism fasta',
                         required=True)
    combine.add_argument('--collection',
                         help='Path to graph collection (text version)',
                         required=True)
    combine.add_argument('--collection-sep',
                         help='separator of the collection file',
                         default='\t')
    combine.add_argument('--collection-selection',
                         help='which graphs to consider in the combination',
                         nargs='+', required=True)
    combine.add_argument('--homology',
                         help='Path to homology graph (text version)',
                         required=True)
    combine.add_argument('--homology-sep',
                         help='separator of the homology file',
                         default='\t')
    combine.add_argument('--seed', help='Path to seed file (text version)',
                         required=True)
    combine.add_argument('--seed-sep', help='separator of the seed file',
                         default='\t')
    combine.add_argument('--seed-threshold',
                         help='defines the threshold of the seed file',
                         default=0.2, type=float)
    combine.add_argument('--output', help='Path to desired output file',
                         required=True)

    diffuse = subparsers.add_parser(
        'diffuse',
        description='S2F diffuse: run the label propagation ' +
                    'on manually provided files',
        help='diffuse command'
    )
    diffuse.set_defaults(func=commands.diffuse)
    diffuse.add_argument('--fasta', help='Path to the organism fasta',
                         required=True)
    diffuse.add_argument('--graph', help='Path to graph (text version)',
                         required=True)
    diffuse.add_argument('--graph-sep', help='separator of the graph file',
                         default='\t')
    diffuse.add_argument('--labelling',
                         help='Path to seed file (text version)',
                         required=True)
    diffuse.add_argument('--labelling-sep', help='separator of the seed file',
                         default='\t')
    diffuse.add_argument('--diffusion-method',
                         help='choose a diffusion method',
                         default='s2f', choices=['s2f', 'consistency-method',
                                                 'label-weighted'])
    diffuse.add_argument('--kernel-parameters', help='kernel parameters',
                         default={'lambda': 0.1}, type=json.loads)
    diffuse.add_argument('--output', help='Path to desired output file',
                         required=True)

    seed_from_hmmer = subparsers.add_parser(
        'hmmer-seed',
        description='S2F hmmer-seed: manually feed ' +
                    'an evalue file and get a seed back',
        help='hmmer-seed command'
    )
    seed_from_hmmer.set_defaults(func=commands.seed_from_hmmer)
    seed_from_hmmer.add_argument('--evalue-file',
                                 help='Path to the evalue file', required=True)
    seed_from_hmmer.add_argument('--obo',
                                 help='Path to obo file (text version)',
                                 required=True)
    seed_from_hmmer.add_argument('--threshold',
                                 help='threshold value to use ' +
                                      'to discard annotations',
                                 default=1e-4, type=float)
    seed_from_hmmer.add_argument('--output',
                                 help='Path to desired output file',
                                 required=True)

    combine_seeds = subparsers.add_parser(
        'combine-seeds',
        description='Aggregates seed files using a linear combination',
        help='combine-seeds command'
    )
    combine_seeds.set_defaults(func=commands.combine_seeds)
    combine_seeds.add_argument('--seed-files', help='seed files to combine',
                               nargs='+', required=True)
    combine_seeds.add_argument('--seed-separator',
                               help='separators of the seed files',
                               default='\t')
    combine_seeds.add_argument('--coefficients',
                               help='list of coefficients, provided in the '
                                    'same order as the seed files, by defaults'
                                    ' it assigns the same value to each seed',
                               nargs='+', default='infer')
    combine_seeds.add_argument('--output', help='Path to desired output file',
                               required=True)

    build_goa_clamp = subparsers.add_parser(
        'build-clamp',
        description='S2F build-clamp: given a fasta file, extracts all the GO annotations from the UniProtKB that'
                    'use the provided evidence codes. This command will fail if run previously to a full S2F '
                    'installation using the `Install` command',
        help='build-clamp command'
    )
    build_goa_clamp.set_defaults(func=commands.build_clamp)
    build_goa_clamp.add_argument('--evidence-codes',
                                 help='manually provide a list of evidence codes that will be used to filter the '
                                      'UniProt GOA. If not provided, S2F will use experimental evidence codes.',
                                 default='experimental')
    build_goa_clamp.add_argument('--fasta', help='Path to the organism fasta', required=True)
    build_goa_clamp.add_argument('--output',
                                 help='Path to desired output file. The fasta filename with an added `.clamp` extension'
                                      'will be created by S2F by default.', default='infer')
    build_goa_clamp.add_argument('--config-file',
                                 help='location of the installation configuration file that will be loaded. If not '
                                      'provided, the default configuration file will be loaded',
                                 default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 's2f.conf'))


    rescore_continuous = subparsers.add_parser(
        'rescore-continuous',
        description='S2F continuous re-scoring: Re-scores the output of S2F to meet the CAFA challenge conditions.',
        help='rescore-continuous command'
    )
    rescore_continuous.set_defaults(func=commands.rescore_continuous)
    rescore_continuous.add_argument('--obo',
                                    help='Path to obo file (text version)',
                                    required=True)
    rescore_continuous.add_argument('--prediction',
                                    help='Path to prediction file (text version)',
                                    required=True)
    rescore_continuous.add_argument('--th-start', type=float,
                                    help='Initial threshold')
    rescore_continuous.add_argument('--outdir', 
                                    help='output_directory')

    args = parser.parse_args()
    if args.subcommand == 'predict':
        if args.run_config is None and (args.alias is None or args.obo is None
                                        or args.fasta is None):
            parser.error('If no run configuration is provided, '
                         'you must provide the following arguments:\n'
                         '\t--alias\n'
                         '\t--obo\n'
                         '\t--fasta')

    try:
        args.func(args)
    except AttributeError as e:
        print(e)
        parser.parse_args(['--help'])
