
# coding: utf-8

# In[1]:


###########################################################################
# Locus Tag Mapper - LTmapper.py
#
# Maps a list of locus tags of different species to your species of choice
# 
# created by Andreas U. Mueller, 2020
###########################################################################
#
#
# What should the program be capable of?
# Batch conversion of locus tags from multiple species to one target species
# What is it NOT for?
# Although LTMapper can do it, it's not intended to look up just single locus tags.
# This can be done easily using existing online tools & databases.
#
#
# TODO
# - user-provided target/query genomes / existing file import
# - generate template config file (with current settings) on the fly if it does not exist (for future runs)
#
# Changelog:
# -- v0.4 - 20201003
#  - timestamped log (to match result files)
#  - revised final statistics output
#  - fixed: lt_getid() accepted partial matches leading to falsely assigned sequences/IDs
# -- v0.3 - 20201001
#  - correct interpretation of no hits BLAST results
#  - fixed query CLI argument
#  - various minor code improvements
# -- v0.2 - 20200926
#  - added CLI parameters
#  - reconciles CLI with config file inputs (CLI takes precedence over config)
#  - reports more statistics
#  - timestamped results
# -- v0.1 - 20200506
#  - initial version, basic functionality
#
###########################################################################

ltm_version = 'v0.4'

import logging
import configparser
config = configparser.ConfigParser()
import os
import sys
import re
import subprocess
import datetime
import pandas as pd


from Bio import SeqIO # pip install biopython
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Bio.Blast import NCBIXML
from io import StringIO


# In[2]:


# make output files unique (e.g. logfile, result tables, ...)
run_timestamp = str(datetime.datetime.now().timestamp()).replace('.', '-')


# In[3]:


# set up CLI arguments
# CLI arguments should override config file settings -> easy adjustments and temporary changes
# Note: defaults must be handled in conjunction with config file parsing (further below),
# otherwise it's not possible to recognize whether an argument has been passed or not (exception: config file name)
import argparse

aparser = argparse.ArgumentParser(
    prog='Locus Tag Mapper '+ltm_version,
    description='Maps locus tags of different species to your species of choice.',
)
# general program behaviour
aparser.add_argument('--version', action='version', version='%(prog)s '+ltm_version, help='display version number')
aparser.add_argument('--overwrite', action='store_true',# default=False,
                     help='overwrite previous files (genomes, BLAST database, CDS FASTA files, ...)')
aparser.add_argument('-v', '--verbose', action='store_true', help='increase verbosity')
#aparser.add_argument('-h', '--help', nargs='?', # this is a default option to argparse and doesn't need to be added
#                     help='print help message') #parser.print_help()

# input control
aparser.add_argument('--conf', nargs='?', default='LTMapper.conf',
                     help='name of configuration file (default: LTMapper.conf)')
aparser.add_argument('-t', '--target', nargs='?',
                     help='target genome (NCBI accession number)')
#aparser.add_argument('-s', '--supplement', nargs='?', #NOT IMPLEMENTED YET
#                     help='supply file with additional sequences for target pool (FASTA format)')
aparser.add_argument('-q', '--query', nargs='*', action='append',
                     help='specify query: -q NCBI_accession_number locus_tag1 locus_tag2 locus_tag2 ...\n'
                     +'OR -q NCBI_accession_number path_to_file\n'
                     +'File must contain one locus tag per line.\n'
                     +'Supply argument multiple times to search for locus tags from different species.')

#aparser.add_argument('-q', '--query_genome', nargs='?', #DEPRECATED
#                     help='search against this genome (NCBI accession number)')
#aparser.add_argument('-l', '--locus_tags', action='append', nargs='+', #DEPRECATED
#                     help='takes multiple values/locus tags OR path to file containing one locus tag per line')

# output control
aparser.add_argument('--workdir', nargs='?',# default='LTM_work/',
                     help='path to work directory, will be created if it does not exist (default: LTM_work/)')
aparser.add_argument('--outputdir', nargs='?',# default='LTM_output/',
                     help='path to output directory, will be created if it does not exist (default: LTM_output/)')
aparser.add_argument('--log', nargs='?', default='LTMapper'+run_timestamp+'.log',
                     help='name of log output file (default: LTMapper_TIMESTAMP.log)')

#aparser.add_argument('-o', '--outfile', nargs='?',# default='LTMapper.out', # there is no single output file...
#                     help='name of output file (default: LTMapper.out)')


# NCBI tools control
aparser.add_argument('--ncbi_email', nargs='?',
                     help='provide your email address to use NCBI web services')
aparser.add_argument('--ncbi_api_key', nargs='?',
                     help='provide your NCBI API key hash to increase NCBI server allowance')

# BLAST control
aparser.add_argument('--path2blast', nargs='?',# default='/usr/bin',
                     help='path to BLAST binaries (default: /usr/bin)')
aparser.add_argument('-e', '--blast_eval', nargs='?',# default='0.0001',
                     help='e-value limit for BLAST search (default: 0.0001)')
aparser.add_argument('--blast_max_target_seqs', nargs='?',# default='1',
                     help='number of BLAST hits to return (default: 1)')

#USE THIS WHEN FINISHED, but currently the work-around below is required
#ltm_ns = aparser.parse_args() # resolve arguments

try: # work-around for invalid arguments; jupyter notebooks pass additional arguments that would be invalid
    # NOTE: if not removed, it will result in having the help message printed twice (stdout and stderr)
    ltm_ns = aparser.parse_args()
except SystemExit:
    ltm_ns = aparser.parse_known_args()[0]
    exc = sys.exc_info()[1]
    print(exc)

#////DEBUG
#aparser.print_help()
#aparser.parse_args(['-v'])
#print(ltm_ns[0].verbosity)
#////DEBUG


# In[4]:


# set up logging
import logging

logfile = ltm_ns.log # from argparse

logger = logging.getLogger('def_logger')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler(logfile, mode='w')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
if ltm_ns.verbose: # change console output based on CLI argument
    ch.setLevel(logging.INFO)
else:
    ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
fh_formatter = logging.Formatter('%(asctime)s (%(levelname)s) %(message)s')
ch_formatter = logging.Formatter('%(levelname)s: %(message)s')
fh.setFormatter(fh_formatter)
ch.setFormatter(ch_formatter)
# add the handlers to logger
logger.addHandler(ch)
logger.addHandler(fh)

# 'application' code
#logger.debug('debug message')
#logger.info('info message')
#logger.warning('warn message')
#logger.error('error message')

#logger.critical('critical message')


# In[5]:


# set up program configuration
# CLI arguments must override configuration file values, handling of no configuration file, missing values, etc

# fetch configuration file name
config_file = ltm_ns.conf # from argparse

# import config
cfread = config.read(config_file)
if cfread:
    logger.info('Importing configuration from file "'+str(config_file)+'" ...\n'
                +'Note: Command line arguments will override configuration file settings.')
else:
    logger.info('Configuration file "'+str(config_file)+'" not found, using command line arguments only ...')

"""#OBSOLETE -> CLI can supplement missing config file values
# verify config file layout
if len(config.sections()) == 0:
    msg = 'Configuration file "'+str(config_file)+'" empty!'
    logger.error(msg)
    sys.exit(msg)
    #raise RuntimeError(msg)
elif not all(sec in config.sections() for sec in ['main', 'NCBI', 'query1']): # verify sections are present
    msg = 'Essential configuration missing! Please adjust the provided template file as intended.'
    logger.error(msg)
    sys.exit(msg)
    #raise RuntimeError(msg)
"""

def resolve_CLI_config(argname, cf_section='', optdefault='', loghandler='', errmsg='', exit=False):
    """
    Resolve CLI argument vs config file argument conflicts. CLI values will override config file.
    CAUTION: Handles only strings.
    When using boolean values, they must be passed as string to optdefault and eval() outside of function.
    """
    flag = eval('ltm_ns.'+argname) # if not supplied by CLI, this will be None
    if not flag:
        try:
            flag = config.get(cf_section, argname)
        except (configparser.NoOptionError, configparser.NoSectionError) as e:
            flag = optdefault # set default value
            if loghandler: loghandler.error(errmsg)
            if exit: sys.exit(errmsg)
    return str(flag)

# overwrite outputs from previous runs
# -> no re-use of genbank/FASTA files, databases, etc; everything will be created from scratch
overwrite_flag = eval(resolve_CLI_config('overwrite', 'main', 'False'))

# set NCBI environment variables
def validEmail(email):
    if len(email) > 7:
        if re.match('^.+@[a-zA-Z0-9-\.]+\.[a-zA-Z]{2,3}$', email) != None:
            return True
        return False

ncbi_email_addr = resolve_CLI_config(
    'ncbi_email', 'NCBI', '', logger,
    'Email address for NCBI Entrez not set! Please provide your email address.',
    True
)

if not validEmail(ncbi_email_addr):
    msg = 'Invalid email address for NCBI Entrez ('+str(ncbi_email_addr)+')! Please provide your email address.'
    logger.error(msg)
    sys.exit(msg)

Entrez.email = ncbi_email_addr

ncbi_key_hash = resolve_CLI_config('ncbi_api_key', 'NCBI')
if ncbi_key_hash:
    Entrez.api_key = ncbi_key_hash
else:
    logger.info('[Optional] NCBI API key hash not set. To increase the number of possible queries per second,'
                    +' provide your key hash (available from your NCBI account).')

# NCBI BLAST tools path
ncbi_btpath = resolve_CLI_config('path2blast', 'main', '/usr/bin')
ncbi_tools = ['makeblastdb', 'blastp'] # essential BLAST programs for LTM

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program): return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file): return exe_file
    return None

# Find BLAST binaries / check if they are in $PATH
# tool version is not checked at the moment! probably requires subprocess launch with -version
if all([which(os.path.join(ncbi_btpath, t)) for t in ncbi_tools]):
    pass
elif all([which(t) for t in ncbi_tools]):
    ncbi_btpath = ''
else:
    msg = 'BLAST exectuables not found. Please make sure to install NCBI BLAST tools correctly.'
    logger.info(msg)
    sys.exit(msg)

# BLAST parameters
b_eval = float(resolve_CLI_config('blast_eval', 'main', '0.0001'))
b_maxts = int(resolve_CLI_config('blast_max_target_seqs', 'main', '1'))

# set up directories
workdir = resolve_CLI_config('workdir', 'main', 'LTM_work/')
outdir = resolve_CLI_config('outputdir', 'main', 'LTM_output/')

# load target species accession from config
target_species = resolve_CLI_config('target', 'main')

# select query settings from config
qry_flag = resolve_CLI_config('query', '', '')
if not qry_flag:
    qry_sec = [qs for qs in config.sections() if qs.startswith('query')]
    qry_bycfg = True
    if not qry_sec:
        msg = 'No query input in configuration file '+str(config_file)+' recognized (check formatting?).'
        logger.error(msg)
        sys.exit(msg)
else:
    qry_sec = eval(qry_flag)
    qry_bycfg = False

# TODO Establish remaining defaults
"""
'-s', '--supplement', nargs='?', #NOT IMPLEMENTED YET
"""

#////DEBUG
#for sec in config.sections():
#    print('section:', sec)
#    for opt in config.options(sec):
#        print('   ', opt, '=', config.get(sec, opt), '---', type(config.get(sec, opt)))
#////DEBUG


# In[6]:


def lt_mkdir(lt_path, loghandler=''):
    """
    Create or re-use directories and log actions.
    """
    try:
        os.makedirs(lt_path.strip('/')) # recursively create given path
        if loghandler: loghandler.info('Creating directory path "'+lt_path+'" ...')
    except FileExistsError:
        # unless the last leaf directory in path already exists
        if loghandler: loghandler.info('Re-using existing directory "'+lt_path+'" ...')
    return None

# set up work directory
#workdir = config.get('main','workdir') #defined above

lt_mkdir(workdir, logger)


# In[7]:


def lt_cds_features(gbrecord, loghandler=''): # current handling of multiple protein IDs produces an intentional error
    """
    Extracts genuine CDS features from genome Biopython GenBank object.
    """
    cds_total = [f for f in gbrecord.features
                 if f.type == 'CDS']
    #genuine_cds = [f for f in gbrecord.features
    #               if f.type == 'CDS' and 'protein_id' in f.qualifiers.keys()]
    genuine_cds = [f for f in cds_total
                   if 'protein_id' in f.qualifiers.keys()]
    # pseudogenes are of type CDS, but don't contain protein_id
    
    # test for ambiguous seq entries
    # NOTE: According to the GenBank flatfile definition, protein_id and translation should be unique to each feature
    aid = [i.qualifiers['protein_id'] for i in genuine_cds if len(i.qualifiers['protein_id']) > 1]
    atr = [tr.qualifiers['protein_id'] for tr in genuine_cds if len(tr.qualifiers['translation']) > 1]
    if len(aid) > 0 or len(atr) > 0:
        errmsg = 'Ambiguous feature qualifiers detected! Check genome '+str(gbrecord.id)+'\n\n'
        +'Ambiguous protein IDs:\n'+str(aid)+'\n\n'
        +'Ambiguous translations (given by protein ID):\n'+str(atr)+'\n'
        
        if loghandler: loghandler.error(errmsg)
        raise RuntimeError(errmsg)
    
    if loghandler: loghandler.info('\nSanity check: FASTA from GenBank for genome '+str(gbrecord.id)
                                   +'\n total features: '+str(len(gbrecord.features))
                                   +'\n total CDS features: '+str(len(cds_total))
                                   +'\n CDS features with protein id: '+str(len(genuine_cds))
                                  )
        
    return genuine_cds, cds_total

def cds2prot(gcdsfeat):
    """
    Convert CDS features to FASTA from genome GenBank record.
    Returns list of Biopython sequence objects (id <-> seq).
    gcdsfeat := genuine CDS features of a GenBank record
    """
    # generate seq records
    # NOTE: protein_id is a unique feature qualifier and contains only one value
    # (but due to parsing by Biopython it's in a list and needs to be retrieved by index 0)
    prot = [SeqRecord(Seq(feat.qualifiers['translation'][0]), id=feat.qualifiers['protein_id'][0])
            for feat in gcdsfeat]
    return prot

def lt_cdsfeat2df(gcdsfeat):
    """
    Returns a pandas dataframe to associate protein id with any possible locus tags.
    """
    fin = []
    
    for feat in gcdsfeat:
        row = [feat.qualifiers['protein_id'][0]]
        
        try:
            row.append(', '.join(feat.qualifiers['locus_tag']))
        except:
            row.append('')
        
        try:
            row.append(', '.join(feat.qualifiers['old_locus_tag']))
        except:
            row.append('')
        
        fin.append(row)
    
    return pd.DataFrame(fin, columns=['protein_id', 'locus_tag', 'old_locus_tag'])

def lt_establish_genome(wdir, acc, loghandler=''):
    """
    wdir := working directory
    acc := NCBI genome accession number
    
    dependencies:
    lt_cds_features
    cds2prot
    lt_cdsfeat2df
    """
    # set up genome work directory
    gdir = os.path.join(wdir, acc)
    lt_mkdir(gdir, logger)
    
    # load genome
    gbkfile = os.path.join(gdir, acc+'.gbk')
    
    if overwrite_flag or not os.path.exists(gbkfile):
        # fetch from NCBI if not existing or overwrite is enforced
        if loghandler: loghandler.info('Fetching genome '+acc+' from NCBI database ...')
        with Entrez.efetch(db='nuccore', id=acc, rettype='gbwithparts', retmode='text') as handle:
            genome = SeqIO.read(handle, 'genbank') # needed for return
        SeqIO.write(genome, gbkfile, 'genbank')
        if loghandler: loghandler.info('Genome file written to '+gbkfile)
    elif os.path.exists(gbkfile):
        # load from disk
        if loghandler: loghandler.info('Using exisiting genome file ('+gbkfile+')')
        genome = SeqIO.read(gbkfile, 'genbank') # needed for return
    #else: # do not catch other error (unless printing the original error traceback for troubleshooting...)
    #    raise RuntimeError('Unexpected error while creating/accessing target genome file '+tgbk)
    
    # validate user-defined against NCBI accession number
    if acc not in genome.id:
        if loghandler: loghandler.warning('Accession number mismatch detected! User-defined genome '
                       +acc+' != '+str(genome.id)+' (NCBI GenBank file)')

    # extract CDS features
    g_cds, g_totalcds = lt_cds_features(genome, loghandler)
    
    # dump fasta using protein ID as sequence headers
    gfastafile = os.path.join(gdir, acc+'_fullcds.fasta')
    
    if overwrite_flag or not os.path.exists(gfastafile):
        gfasta = cds2prot(g_cds) # convert
        nseqdump = SeqIO.write(gfasta, gfastafile, 'fasta') # dump to disk
        if loghandler: loghandler.info('Extracted '+str(nseqdump)+' sequences and wrote to '+gfastafile)
    elif os.path.exists(gfastafile):
        gfasta = [s for s in SeqIO.parse(gfastafile, 'fasta')] # needed for return
        gcds = lt_cds_features(genome)[0]
        if loghandler: loghandler.info('Loaded '+str(len(gfasta))+' sequences from existing file '+gfastafile)
    #else: # do not catch other error (unless printing the original error traceback for troubleshooting...)
    #    raise RuntimeError('Unexpected error while creating/accessing '+tfasta)
    
    # create dataframe for protein id <-> locus tag reference
    g_lt_ref = lt_cdsfeat2df(g_cds)
    
    return genome, gbkfile, gfasta, gfastafile, g_lt_ref


# In[8]:


### prepare target species

# load target species accession from config
#target_species = config.get('main', 'target_species') # handled above

# establish target genome basic data
logger.info('Setting up target genome data ...')
tgenome, tgbk, tfasta, tfastafile, target_lt_ref = lt_establish_genome(workdir, target_species, logger)

'''
# add supplement sequences
try:
    tsfile = config.get('main','target_supplement')
except configparser.NoOptionError:
    logger.info('No target supplement specified, continuing ...')
    pass

def lt_add_supplement_seq(suppl_file, ):
    return None

if os.path.exists(tsfile):
    # load supplemental sequences into list object
    tsuppl = [s for s in SeqIO.parse(tsfile, 'fasta')]
    logger.info('Target supplement specified. Loaded '+str(len(tsuppl))+' additional sequences from '+tsfile)
    # resolve conflicts; mode: override existing entries by id
    fl_pre = [sr for sr in tfasta if sr.id not in [i.id for i in tsuppl]]
    fl = fl_pre+tsuppl
    tfasta_fin = os.path.join(os.path.split(tgbk)[0], target_species+'_cds+suppl.fasta')
    nseqdump = SeqIO.write(fl, tfasta_fin, 'fasta') # dump to disk
    logger.info('Wrote '+str(nseqdump)+' sequences to '+tfasta_fin
                +'\nFrom genome: '+str(len(tfasta))
                +'\nSupplemented: '+str(len(tsuppl))
                +'\nOverridden by supplement: '+str(len(tfasta)-len(fl_pre))
               )
elif not os.path.exists(tsfile):
    tfasta_fin = tfasta
    logger.info('No target supplement found or specified, continuing ...')
else:
    raise RuntimeError('Unexpected error while accessing '+tsfile)
'''

tfasta_fin = tfasta # how to supplement sequences? is it necessary?
tfastafile_fin = tfastafile

if not len(tfasta_fin) > 0:
    logger.error('Unexpected error. No sequences in FASTA file of target genome!')
    sys.exit(1)

# generate local database
localblastdb = os.path.join(os.path.split(tgbk)[0], target_species)
if not all([os.path.exists(dbfile) for dbfile
            in [localblastdb+'.psq', localblastdb+'.phr', localblastdb+'.pin']]):
    logger.info('Creating local BLAST database for target species '+target_species+' from '+tfastafile_fin+' ...')
    mkdb_cline = os.path.join(ncbi_btpath, 'makeblastdb')+' -dbtype prot -in '+tfastafile_fin+' -title '+target_species+' -out '+localblastdb
    mkdb_proc = subprocess.Popen(mkdb_cline.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #set PIPE to DEVNULL to discard
    mkdb_output, mkdb_error = mkdb_proc.communicate()
    # report outcome to log
    if mkdb_error.decode('utf-8') != '':
        logger.error('Errors occurred during makeblastdb runtime:\n'+mkdb_error.decode('utf-8'))
    if mkdb_output.decode('utf-8') != '':
        logger.info('makeblastdb output:\n'+mkdb_output.decode('utf-8'))
else:
    logger.info('Using existing local BLAST database in '+localblastdb)


# In[9]:


def lt_getid(query_locus_tags, refdf):
    """
    Get protein id for query locus tags.
    """
    fq_mapping = {}
    refmap = list(zip(refdf['protein_id'], refdf['locus_tag'], refdf['old_locus_tag']))
    
    # CAVEATS: case-sensitive, overwrites value in case of duplicates
    for q in query_locus_tags:
        fq_mapping.update(
            #[f(x, y) for x, y in zip(df['col1'], df['col2'])]
            #{q:row[0] for row in refmap if q in row[1] or q in row[2]}
            # PROBLEM DISCOVERED: this will accept half-matches (ie. Rv0063 matches Rv0063a)
            # attempted solution:            
            {q:row[0] for row in refmap
             if any([i == q for i in row[1].split(', ')])
             or any([i == q for i in row[2].split(', ')])}
        )
    
    return fq_mapping


# In[10]:


## prepare queries

# select query settings from config
#qry_sec = [qs for qs in config.sections() if qs.startswith('query')] # handled above

def fetch_queryLT_from_file(qfile):
    # load user-given locus tags
    with open(qfile, 'r') as handle:
        qlts = [l.strip() for l in handle.readlines() if l.strip()]
    return qlts

# proxy to collect results for summary
result_collector = []

# initialize counter for statistics
# this will become a nested dictionary; top level for query, sub-levels for different the countings
counter = {}

# handle each query set
for query in qry_sec:
    # set up counter for statistics
    counter[query] = {}
    counter[query]['mapped'] = 0
    counter[query]['unmapped'] = 0
    counter[query]['noid'] = 0
    counter[query]['noseq'] = 0
    counter[query]['noblast'] = 0
    
    # determine input method for query
    if qry_bycfg: # variable from CLI/config resolution code above
        qacc = config.get(query, 'query_genome')
        qlt_file = config.get(query, 'locus_tags')
        qlts = fetch_queryLT_from_file(qlt_file)
    else:
        qacc = query[0]
        qlt_file = query[1]
        if os.path.exists(qlt_file):
            qlts = fetch_queryLT_from_file(qlt_file)
        else:
            qlts = query[1:]
    
    logger.info('Query: received '+str(len(qlts))+' locus tags for genome '+str(qacc))
    
    # establish query genomes
    qgenome, qgbk, qfasta, qfastafile, query_lt_ref = lt_establish_genome(workdir, qacc, logger)
    
    qdir = os.path.split(qgbk)[0] # set query directory
    
    # associate locus tags with protein id (ie handle to access seq later)
    qid_file = os.path.join(qdir, qacc+'_lt2id.tsv')
    
    #if not os.path.exists(qid_file) or overwrite_flag: # is this saving so much time?
    # query selection will be different in most cases and thus should always be freshly written
    qid_map = lt_getid(qlts, query_lt_ref)
    qid_leftover = [v for v in qlts if not v in qid_map.keys()]
    
    counter[query]['noid'] = len(qid_leftover)
    
    logger.info('Locus tag <-> protein ID: '+str(len(qid_map))+'/'+str(len(qlts)))
    
    # write association file to disk
    # query selection will be different in most cases and thus should always be freshly written
    # does it need to be written?
    if qid_map:
        qid_map_df = pd.DataFrame.from_dict(qid_map, orient='index') # data frame conversion only because of file output
        qid_map_df.reset_index(inplace=True)
        qid_map_df.rename(columns={'index':'locus_tag', 0:'protein_id'}, inplace=True)
        qid_map_df.to_csv(qid_file, sep='\t', index=False)
    else:
        logger.warning('None of the locus tags were found in genome '+str(qacc)+'!')
        continue
    
    if qid_leftover: # inform user
        qid_leftover_msg = ''
        for i in qid_leftover:
            qid_leftover_msg += i+'\n'
        logger.warning('Locus tags without ID:\n'+str(qid_leftover_msg))
    
    if len(qid_map) != len(qid_map_df):
        logger.warning('Conversion of locus tag <-> protein id association to data frame created mismatch!'
                       +'\n#before '+str(len(qid_map))+' #after '+str(len(qid_map_df))
                      )
    
    #elif os.path.exists(qid_file):
    #    logger.info('Found '+str(qid_file)+' ... loading from disk ...')
    #    qid_map_df = pd.read_csv(qid_file, sep='\t').fillna('')
    #else: # do not catch other error (unless printing the original error traceback for troubleshooting...)
    #    raise RuntimeError('Unexpected error while accessing '+qid_file)
    
    # select sequences
    #qselect = [sr for sr in qfasta if sr.id in qid_map.values()]
    # for some reason, this does not yield all valid records
    # update: this was most likely related to a wrong matching in lt_getid()
    # current workaround
    qselect = []
    qselect_leftover = []
    qfasta_bykey = {sr.id:sr for sr in qfasta}
    for sel in qid_map:
        try:
            qselect += [qfasta_bykey[qid_map[sel]]]
        except KeyError:
            qselect_leftover += [(sel, qid_map[sel])]
    #qselect_leftover = [(i, qid_map[i]) for i in qid_map if qid_map[i] not in [sr.id for sr in qselect]]
    
    counter[query]['noseq'] = len(qselect_leftover)
    
    logger.info('Locus tag/ID <-> sequence: '+str(len(qselect))+'/'+str(len(qid_map)))
    
    if not qselect:
        logger.error('Unexpected error: Assignment of locus tags to sequences returned 0 sequences.')
        continue
    
    if qselect_leftover: # inform user
        qselect_leftover_msg = ''
        for i in qselect_leftover:
            qselect_leftover_msg += i+'\n'
        logger.warning('Locus tags/IDs without sequence:\n'+str(qselect_leftover_msg))
    
    # write selected seq to file for BLAST
    # actually not needed, can run BLAST from STDIN
    # in any case it has to be written every time, as each query is different is most cases
    qfasta_sel = os.path.join(qdir, qacc+'_cds_query.fasta')
    
    #if overwrite_flag or not os.path.exists(qfasta_sel):
    nseqdump = SeqIO.write(qselect, qfasta_sel, 'fasta') # dump to disk
    logger.info('Selected '+str(nseqdump)+' sequences for query and wrote to '+qfasta_sel)
    #elif os.path.exists(qfasta_sel):
    #    qselect = list(SeqIO.parse(qfasta_sel, 'fasta')) # load from disk (needed for BLAST)
    #else: # do not catch other error (unless printing the original error traceback for troubleshooting...)
    #    raise RuntimeError('Unexpected error while creating/accessing '+qfasta_sel)
    
    # blast against target (local db)
    #b_eval = config.getfloat('main', 'blast_eval') # defined above
    #b_maxts = config.getint('main', 'blast_max_target_seqs') # defined above
    # localblastdb is defined when makeblastdb is run (see above)
    
    logger.info('Mapping '+str(nseqdump)+' locus tags from genome '+str(qgenome.id)+' to target '+str(tgenome.id))
    
    #construct BLAST command
    blastp_cline = os.path.join(ncbi_btpath, 'blastp')+' -db '+localblastdb+' -outfmt 5 -evalue '+str(b_eval)+' -max_target_seqs '+str(b_maxts)
    
    qresult = pd.DataFrame() # create dataframe proxy to collect BLAST results
    
    # loop over query sequences and search each against target genome using BLAST
    for qseqrec in qselect:
        # define command-line process
        blastp_proc = subprocess.Popen(
            blastp_cline.split(),
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        
        # transmit query sequence to BLAST process
        blastp_out, blastp_err = blastp_proc.communicate(str(qseqrec.seq).encode('utf-8'))
        
        blastp_out = blastp_out.decode('utf-8')
        blastp_err = blastp_err.decode('utf-8')
        
        # report errors to log
        if blastp_err:
            logger.error('Errors occurred during blastp run of '+q+':\n'+blastp_err)
        
        # parse BLAST result
        blastp_result = NCBIXML.read(StringIO(blastp_out))
        
        if blastp_result.descriptions:
            # construct result dataframe
            blastp_df = pd.DataFrame([
                [rec.title.split(' ')[1], rec.score, rec.e]
                for rec in blastp_result.descriptions], columns=['target_id', 'score','eval'])
            counter[query]['mapped'] += 1
        else:
            # handle no hits result here
            blastp_df = pd.DataFrame([['', pd.np.nan, pd.np.nan]], columns=['target_id', 'score','eval'])
            counter[query]['noblast'] += 1

        # add query id to be able to merge
        blastp_df['protein_id'] = qseqrec.id
        
        # merge results with final table
        qresult = qresult.append(blastp_df)
    
    # compile output result for current query
    qresult = qresult.merge(query_lt_ref, on='protein_id', how='left') # add query locus tags
    qresult = qresult.merge(
        target_lt_ref.rename(
            columns={'protein_id':'target_id',
                     'locus_tag':'target_locus_tag',
                     'old_locus_tag':'target_old_locus_tag'}
        ), on='target_id', how='left') # add target locus tags
    
    qresult = qresult[
        ['protein_id', 'locus_tag', 'old_locus_tag',
         'target_id', 'target_locus_tag', 'target_old_locus_tag', 'score', 'eval']] # sort columns
    
    # write final result to disk
    lt_mkdir(outdir) # will try to create or re-use exisiting directory
    
    qfinalfile = os.path.join(outdir, qacc+'_result_'+run_timestamp+'.tsv')
    
    qresult.to_csv(qfinalfile, sep='\t', index=False)
    
    # summary statistics for current query results
    # will only work for max_target_seqs = 1, otherwise it becomes very complex -> should use counter during iteration
    #num_mapped = len(qresult[qresult['target_id'].astype(bool)].drop_duplicates['protein_id'])
    #num_duplicated = len(qresult[qresult['protein_id'].duplicated()])
    #num_unmapped = len(qresult[qresult['target_id'].astype(bool)!=True])
    counter[query]['unmapped'] = counter[query]['noblast'] + counter[query]['noid'] + counter[query]['noseq']
    
    logger.info(
        '\nSummary: '
        +str(qgenome.id)+' ('+str(qgenome.description.split(',')[0])[:50]+') -> target '
        +str(tgenome.id)+' ('+str(tgenome.description.split(',')[0])[:50]+')'
        +'\n Successful mappings: '
        +str(counter[query]['mapped'])+' of '+str(len(qlts))+' (total input)'
        +'\n Total result entries: '+str(len(qresult))
        #+'\n Ambiguous matches in target: '+str(num_duplicated)
        +'\n Query sequences without BLAST result: '+str(counter[query]['noblast'])
        +'\n Unrecognized locus tags (no ID or sequence): '+str(counter[query]['noid'] + counter[query]['noseq'])
        +'\n Total query locus tags unassigned: '+str(counter[query]['unmapped']) #str(num_unmapped+len(qlts)-len(qselect))
    )
    
    result_collector += [qfinalfile] # not in use yet


# In[11]:


# clean-up and exit
sys.exit(0)


# In[ ]:


### DEBUG
#myrec = SeqRecord(
#    Seq('MAVVPLGEVRNRLSEYVAEVELTHERITITRHGHPAAVLISADDLASIEETLEVLRTPGASEAIREGLADVAAGRFVSNDEIRNRYTAR'), id='myrec')
myrec = SeqRecord(
    Seq('MRSTVAVAVAAAVIAASSGCGSDQPAHKASQSMITPTTQIAGAGVLGNDRKPDESCARAAAAADPGPPTRPAHNAAGVSPEMVQVPAEAQRIVVLSGDQLDALCALGLQSRIVAAALPNSSSSQPSYLGTTVHDLPGVGTRSAPDLRAIAAAHPDLILGSQGLTPQLYPQLAAIAPTVFTAAPGADWENNLRGVGAATARIAAVDALITGFAEHATQVGTKHDATHFQASIVQLTANTMRVYGANNFPASVLSAVGVDRPPSQRFTDKAYIEIGTTAADLAKSPDFSAADADIVYLSCASEAAAERAAVILDSDPWRKLSANRDNRVFVVNDQVWQTGEGMVAARGIVDDLRWVDAPIN'),
    id='myrec')

blastp_proc = subprocess.Popen(
    blastp_cline.split(),
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
)

blastp_out, blastp_err = blastp_proc.communicate(str(myrec.seq).encode('utf-8'))

blastp_out = blastp_out.decode('utf-8')
blastp_err = blastp_err.decode('utf-8')


blastp_result = NCBIXML.read(StringIO(blastp_out))
### DEBUG

