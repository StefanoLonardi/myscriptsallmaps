#!/usr/bin/pypy
#
# pacbio_lg_assignment
#
# This runs a blast of the /blastdb/Cowpea_iSelect_Ref_seq_raw.fa against the pacbio assembled contigs and unassembled.
#
# It then tallys the hits over the e-score threshold to each LG for each contig.

import os, sys, sqlite3

def runit(cmd):
  print cmd
  os.system(cmd)

exp_cutoff = 50   # 1e-50
consensus_map = '/home/stelo/ALLMAPS_in_progress/CI/CI_map_AllSNPs.txt'
contig_stats_dict = {}  # contig_stats_dict[contig] = [seqlen, is_contig]  Will be filled by get_contig_stats()
snp_wgs_dict = {}       # Xref from snp to wgs_snp
consensus_map_dict = {} # given contig, returns [lg, cm]

# vars for sqlite database
conn = 0
c_snps_contigs = 0
c_lg_assign = 0

# Routine to Get counts of bases in each contig, and GC.  (Not Writing to outfile)
# Loads contig_stats_dict[contig] = [contig_bp, is_contig n(1)]
def get_contig_stats(filename, outfile):
  lines = open(filename).read().splitlines()   # 'cowpea.contigs.fasta'  Read the entire file into a list. Trims line endings.

  if 'CONTIG' in filename.upper():
    is_contig = 1
  else:
    is_contig = 0
  num_lines = len(lines)
  print filename, ' Lines: ', num_lines

  #outh = open(outfile, 'w')  # 'cowpea_pacbio_contig_stats.txt'
  #outh.write('#Node\tBP\tNum N\tNon_N_BP\tGC\n')

  startline = 0
  stopline = 0

  statsloaded = 0
  lnum = 0
  while lnum < num_lines:
    header = lines[lnum]
    if header[0] != '>':
      print 'Error. line ', lnum, ' is not a header. ', header
      sys.exit(1)

    contig = header.split(' ')[0][1:]   # >tig00000019 len=34616 read...  --> tig00000019

    # Find the sequence lines that make up this contig
    lnum += 1
    startline = lnum
    while lnum < num_lines:
            if lines[lnum][0] == '>':
                    stopline = lnum
                    break
            lnum += 1
    # The last contig
    if lnum >= num_lines:
      stopline = lnum

    # Make a byte array out of that section of the file.  Can join. (No line endings)
    seq = ''.join(lines[startline:stopline])
    seqlen = len(seq)
    num_n = seq.count('N')
    non_n = seqlen - num_n

    gc_count = seq.count('G') + seq.count('C')
    gcnum = (gc_count * 1.000) / (non_n * 1.000)
    gc = "%.3f" % gcnum  # String with 3 decimals

    contig_stats_dict[contig] = [seqlen, is_contig]
    blank = ''

    c_lg_assign.execute('insert into lg_assign \
     (contig, is_contig, contig_bp, lg1, lg2, lg3, lg4, lg5, lg6, lg7, lg8, lg9, lg10, lg11, snp_tot, top1, top2, pct1, num_to_lg, contig_lg75, contig_lg100, assign75, assign100) \
     values \
     (     ?,         ?,       ?,   0,   0,   0,   0,   0,   0,   0,   0,   0,    0,    0,       0,    0,    0,  0.0,         ?,           ?,            ?,        0,         0)', \
      (contig, is_contig,  seqlen, blank, blank, blank)) 
    conn.commit()
    statsloaded += 1

    #outh.write(contig + '\t' + str(seqlen) + '\t' + str(num_n) + '\t' + str(non_n) + '\t' + gc + '\n')
    if lnum >= num_lines:
      break

  conn.commit()
  #outh.close()
  #print 'Done creating ', outfile
  return
# End of get_contig_stats

# Load a Blast file
def load_blast_file(blastfile):
  print 'Loading ', blastfile
  numloaded = 0
  lnum = 0
  for line in open(blastfile):
    line = line.strip()
    if 'qacc' in line:
      continue  # Throw away the header
    cols = line.split('\t')
    snp         = cols[0].ljust(7)
    contig      = cols[1]  #.ljust(20)
    escore      = cols[2]
    escore      = escore.replace('.00E', 'E')
    exponent    = int(cols[3])
    hitlen      = int(cols[4])
    nident      = int(cols[5])
    frames      = cols[6]
    snp_start   = int(cols[7])
    snp_end     = int(cols[8])
    contig_start = int(cols[9])
    contig_end   = int(cols[10])

    # -- GET PACBIO CONTIG LENGTH, and is_contig (0/1) -- 
    contig_len, is_contig  =  contig_stats_dict[contig]

    # Get wgs_snp
    wgs_snp = ''
    if snp in snp_wgs_dict:
      wgs_snp = snp_wgs_dict[snp]

    # Get map info
    lg_2015 = ''
    cm_2015 = ''
    if snp in consensus_map_dict:
      lg_2015, cm_2015 = consensus_map_dict[snp]

    lnum += 1
   
    if exponent >= exp_cutoff:
      blank = ''
      c_snps_contigs.execute('insert into snps_contigs \
       (snp, wgs_snp, contig, is_contig, contig_len, escore, exponent, hitlen, nident, snp_start, snp_end, contig_start, contig_end, lg_2015, cm_2015, cm_gap, pct1, contig_lg75, contig_lg100, assign75, assign100, min_cm, max_cm, cm_range) \
       values \
       (  ?,       ?,      ?,         ?,          ?,      ?,        ?,      ?,     ?,          ?,       ?,            ?,          ?,      ?,       ?,       ?,    ?,         ?,          ?,        ?,         ?,      ?,     ?,         ?)', \
       (snp, wgs_snp, contig, is_contig, contig_len, escore, exponent, hitlen, nident, snp_start, snp_end, contig_start, contig_end, lg_2015, cm_2015,    0.0,  0.0,     blank,      blank,        0,         0,    0.0,    0.0,      0.0)) 
      conn.commit()
      numloaded += 1

  print blastfile, '  Loaded ', numloaded
  conn.commit()
  return     
# End of load_blast_file

# ------------------------------ START MAIN ROUTIINE --------------------------

if len(sys.argv) < 2:
  input_files = ['cowpea.contigs.fasta', 'cowpea.unassembled.fasta']
  print 'No parameters passed.  Assuming', input_files[0], input_files[1]
else:
  # Make a list of the input files passed.  These are usually the pacbio contigs file and unassembled file.
  input_files = [sys.argv[1]]
  if len(sys.argv) > 2:
    input_files.append(sys.argv[2])  # They want to do the unassembled file also

refseqfile = 'Cowpea_iSelect_Ref_seq_raw.fa'

# Check the input files so we don't have to mess wit the file names as much
for infile in input_files:
  if '/' in infile:
    print 'Please run this in the folder where the pacbio contigs file is.'
    sys.exit(1)
  if not os.path.exists(infile):
    print 'Error.  Could not find input file: ', infile
    sys.exit(1)

# Load snp_wgs_dict to allow looking up WGS_SNP given a SNP
# >1_0002  35064055_01450
# TAAAGTTGAATTTGGCCTTCCTGATTTAGAG
snp_wgs_dict = {}
for line in open('/blastdb/' + refseqfile):
  if line[0] == '>':
    line = line.strip()
    cols = line[1:].split()  # ['1_0002', '35064055_01450']
    snp_wgs_dict[cols[0].ljust(7)] = cols[1]

# Load consensus_map_dict so we can look up lg/cm
# 2_32511 6       34.67
# 2_50677 6       34.67
# 2_33821 6       34.67  
consensus_map_dict = {}
for line in open(consensus_map):
  cols = line.strip().split('\t')   # ['2_32511', '6', '34.67']
  consensus_map_dict[cols[0].ljust(7)] = cols[1:]  # list of lg and cm   Have to pad to 7 spaces for "1_0002 "

# Create sqlite tables

# SNP/CONTIG DETAIL TABLE
conn = sqlite3.connect(':memory:')

# -------------------------
c_snps_contigs = conn.cursor()
c_snps_contigs.execute('CREATE TABLE snps_contigs (snp text, wgs_snp text, contig text, is_contig integer, contig_len integer, escore text, exponent integer, hitlen integer, nident integer, snp_start integer, \
snp_end integer, contig_start integer, contig_end integer, lg_2015 text, cm_2015 text, cm_gap real, pct1 real, contig_lg75 text, contig_lg100 text, \
assign75 integer, assign100 integer, min_cm real, max_cm real, cm_range real)')
snps_contigs_table_info = c_snps_contigs.execute("PRAGMA table_info(snps_contigs)").fetchall()  # (0, u'contig', u'text', 0, None, 0)
snps_contigs_headers  = list(map(lambda x: x[1], snps_contigs_table_info))
snps_contigs_coltypes = list(map(lambda x: x[2], snps_contigs_table_info))

# ASSIGNMENT TABLE (WILL FILL WHEN GETTING CONTIG STATS)
c_lg_assign = conn.cursor()
c_lg_assign.execute('create table lg_assign \
(contig text, is_contig integer, contig_bp integer, lg1 integer, lg2 integer, lg3 integer, lg4 integer, lg5 integer, lg6 integer, lg7 integer, lg8 integer, lg9 integer, lg10 integer, lg11 integer, snp_tot integer, top1 integer, top2 integer, pct1 real, num_to_lg text, contig_lg75 text, contig_lg100 text, assign75 integer, assign100 integer)')
lg_assign_table_info = c_lg_assign.execute("PRAGMA table_info(lg_assign)").fetchall()
lg_assign_headers  = list(map(lambda x: x[1], lg_assign_table_info))
lg_assign_coltypes = list(map(lambda x: x[2], lg_assign_table_info))
  
# Do the Blasts and generate contig stats files
contig_stat_dict = {}    # Will be contig_stat_dict[contig] = [contig_bp, is_contig True/False]

for infile in input_files:
  if os.path.exists('haveblasts.flg'):
    print 'Found haveblasts.flg. Not re-running the Blasts.'
  else:
    runit('formatdb -i ' + infile + ' -p F -o F')
    runit('runBlastn_all.py /blastdb/' + refseqfile + ' ' + infile)
  blastoutfile = refseqfile + '.' + infile 
  contigstatsfile = infile + '_contig_stats.txt'
  get_contig_stats(infile, contigstatsfile)
  load_blast_file(blastoutfile)

# Index lg_assign for update speed
conn.execute('CREATE INDEX lg_assign1 ON lg_assign(contig)')

# Get a summary of the number of hits by contig.  Use this to populate lg_assignment
conn.execute('CREATE INDEX snps_contigs1 ON snps_contigs(contig, lg_2015)')

print 'Sumarize 1...'
conn.execute("create table hitsummary as select contig, lg_2015, count(*) as numhits from snps_contigs where lg_2015 <> '' group by contig, lg_2015")
conn.execute("create index hitsum1 on hitsummary(contig, lg_2015, numhits)")
print 'Sumarize 2...'
c_contigs = conn.cursor()
c_contigs.execute("select distinct contig from hitsummary")
conn.commit()
c_contigs_rows = c_contigs.fetchall()
numContigRows = len(c_contigs_rows)
print 'Contigs with a BLAST hit', numContigRows
#pctHeader =  '---10---20---30---40---50---60---70---80---90--100'
#print pctHeader

c_hitspercontig = conn.cursor()

lpct = 0
c_contigs_rownum = -1
for row in c_contigs_rows:
  c_contigs_rownum += 1
  contig = row[0]
  #pctdone = int((c_contigs_rownum * 1.0) / (numContigRows * 1.0) * 100)  
  #if pctdone > lpct+1:
  #  sys.stdout.write('.')
  #  sys.stdout.flush()
  #  lpct = pctdone

  # Count hits to each LG for a contig
  lg_hits = {}
  for lg in range(1,12):  # 1-11  Create dict and zero it.
    lg_hits[lg] = 0

  #c_hitspercontig.execute("select lg_2015, count(*) as numhits from snps_contigs where contig = ? and lg_2015 <> ? group by lg_2015", (contig, ''))
  c_hitspercontig.execute("select lg_2015, numhits from hitsummary where contig = ?", (contig,))
      
  for lg_2015, num in c_hitspercontig:
    lgnum = int(lg_2015)
    if lgnum in lg_hits:
      lg_hits[lgnum] += num

  # -- GET TOP HIT LG -- 
  mtopLGcount = 0
  mtopLGnum   = 0
  msnp_tot    = 0
  for a in range(1,12):
    if lg_hits[a] > mtopLGcount:
      mtopLGcount = lg_hits[a]
      mtopLGnum = a
    msnp_tot = msnp_tot + lg_hits[a]
  
  # -- GET 2ND HIT LG -- 
  m2ndLGcount = 0
  m2ndLGnum = 0
  for a in range(1,12):
    if a == mtopLGnum:
      continue   # already found the top one.  Skip it. 
    if lg_hits[a] > m2ndLGcount:
      m2ndLGcount = lg_hits[a]
      m2ndLGnum = a

  if msnp_tot == 0:
    mpct = 0
  else:
    mpct = round(mtopLGcount * 100.0 / ((mtopLGcount + m2ndLGcount) * 1.0),1)
  mnum_to_lg = str(mtopLGcount) + '/' + str(msnp_tot)  # 8/9
  mcontig_lg75  = 0
  mcontig_lg100 = 0
  massign75   = 0
  massign100  = 0

  if mpct >= 75:
    mcontig_lg75 = mtopLGnum
    massign75 = 1

    if mpct == 100:
      mcontig_lg100 = mtopLGnum
      massign100 = 1
  
  # -- UPDATE RECORD FOR THIS CONTIG TO LG_ASSIGN
  #                                          
  c_lg_assign.execute('update lg_assign set lg1 = ?, lg2 = ?, lg3 = ?, lg4 = ?, lg5 = ?, lg6 = ?, lg7 = ?, lg8 = ?, lg9 = ?, lg10 = ?, lg11 = ?, \
   snp_tot = ?, top1 = ?, top2 = ?, pct1 = ?, num_to_lg = ?, contig_lg75 = ?, contig_lg100 = ?, assign75 = ?, assign100 = ? \
   where contig=?', (lg_hits[1], lg_hits[2], lg_hits[3], lg_hits[4], lg_hits[5], lg_hits[6], lg_hits[7], lg_hits[8], lg_hits[9], lg_hits[10], lg_hits[11], \
    msnp_tot, mtopLGcount, m2ndLGcount, mpct,    mnum_to_lg, str(mcontig_lg75),  str(mcontig_lg100), massign75, massign100, contig)) 
  conn.commit()
  
  # UPDATE snps_contigs WITH CONTIG ASSIGNMENTS TO LG (created above)
  c_snps_contigs.execute('update snps_contigs \
   set pct1=?, contig_lg75=?, assign75=?, contig_lg100=?,assign100=? where contig=?', \
        (mpct, str(mcontig_lg75), massign75, str(mcontig_lg100), massign100, contig))
  conn.commit()

# DONE WITH LOOP
conn.commit()

# MARK CONTIG MIN AND MAX CM FOR ASSIGN100 SNPs
print 'Mark min/max cM'
c_temp = conn.cursor()
c_temp.execute("SELECT contig, MIN(cast(cm_2015 as real)) as min_cm, MAX(cast(cm_2015 as real)) as max_cm FROM snps_contigs \
  WHERE assign100 = 1 AND cm_2015<>'' GROUP BY contig")

lpct = 0
rownum = 0
for row in c_temp:
  contig = row[0]
  min_cm = row[1]
  max_cm = row[2]
  cm_range = abs(max_cm-min_cm)
  c_snps_contigs.execute('update snps_contigs set min_cm = ?, max_cm = ?, cm_range = ? where contig = ? and assign100 = 1', (min_cm, max_cm, cm_range, contig))
    
# ------- DONE CALCULATING ---------
  
# Output snps_contigs
print 'Output snps_contigs.txt'
mh = open(input_files[0]+'_CI_snps_contigs.csv', 'w')
#mh.write('\t'.join(snps_contigs_headers) + '\n')

#c_snps_contigs.execute('select * from snps_contigs order by snp, contig, contig_start')
c_snps_contigs.execute('select contig, contig_start, lg_2015, cm_2015 from snps_contigs where lg_2015<>"" order by lg_2015, cm_2015, contig_start') 
allrows = c_snps_contigs.fetchall()
print 'snps_contigs rows: ', len(allrows)
rownum = 0
for row in allrows:
  outset = []
  for colnum in range(0, 4): #len(snps_contigs_coltypes)):
    if snps_contigs_coltypes[colnum] == 'integer':
      outset.append(str(row[colnum]))
    elif snps_contigs_coltypes[colnum] == 'real':
      outset.append(str(round(row[colnum],3)))
    else:
      outset.append(str(row[colnum]))
  line = ','.join(outset) + '\n'
  mh.write(line)
  rownum += 1
mh.close()

#print 'Output lg_assign'
####0 (0, u'contig', u'text', 0, None, 0)
####1 (1, u'is_contig', u'integer', 0, None, 0)

#c_lg_assign.execute('select * from lg_assign')
#mh = open('lg_assign.txt', 'w')
#mh.write('\t'.join(lg_assign_headers) + '\n')

#for row in c_lg_assign:
#  outset = []
#  for colnum in range(0, len(lg_assign_coltypes)):
#    if lg_assign_coltypes[colnum] == 'integer':
#      outset.append(str(row[colnum]))
#    elif lg_assign_coltypes[colnum] == 'real':
#      outset.append(str(round(row[colnum],3)))
#    else:
#      outset.append(row[colnum])
#  line = '\t'.join(outset) + '\n'
#  mh.write(line)
#mh.close()
print 'Done.\n'
#os.system('ls -l')
