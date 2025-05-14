#!/usr/bin/python

import math
import string
from numpy import *
from numpy.dual import *
import struct

DEBUG = 0

HDR_SIZE = 1024		#size of file header


class PF:
  p = {} # dict of parameter key and value
  s = [] # ordered list of keys (dict is not keeping sequence)

  def __init__(self, pf):
    # initialize class -> read in parameters from an open file
    p = {}
    s = []
    kv = []
    pv = []
    line = pf.readline(80)	#read up to 80 bytes/line

    while (line):
#      print 'read line: ',line, len(line), ord(line[-1])
      if len(line) == 80 and line[-1] != '\n':	
        break

      if line[0] == '#' or len(line.strip()) == 0:
        line = pf.readline(80)	#read next line
        continue
   
      kv = line.split(":")
      if len(kv) < 2:		#if not a comment or blank line then there must be a :
        print "ERROR: parameter file missing : seperator :",line
        break

      kv[0] = kv[0].strip()
      kv[1] = kv[1].strip()
      p[kv[0]] = kv[1]
      s.append(kv[0])
      pv = kv[1].split(",")

      for i in range(len(pv)):
        pv[i] = pv[i].strip()

      p[kv[0]] = pv
      line = pf.readline(80)	#read next line

    self.p = p # store dict
    self.s = s # store keyword list

  def get(self,key):
    # get value for key
    if self.p.has_key(key):
      return self.p[key]
    else:
      return "" # if key does not exist

  def set(self,key,value):
    # set value for key
    if  self.p.has_key(key) == 0:
      self.s.append(key) # append new key to list 
    self.p[key]=value # write key/value
    
  def dump(self):
    # print all key/value pairs
    for key in self.s:
      print "%s:\t\t %s" % (key, self.p[key])

  def write(self,pf):
    # write key/value pairs separated by : <tab> to open file, do not close file
    tstr = ""
    for key in self.s:
      if (type(self.p[key])) == list:
        for i in range(len(self.p[key])):
          if i == 0:
            tstr = self.p[key][i]
          else :
            tstr = "%s, %s" % (tstr, self.p[key][i])
        pf.write("%s:\t\t %s\n" % (key, tstr))        
      else :
        pf.write("%s:\t\t %s\n" % (key, self.p[key]))        


## (fs, ns, psf, ncycle, cycle_delay, elev0, elev_step, nelev, pmodes, hdr_len, rec_len) = Read_FHeader(df)
#  Read raw data file header from file
#
#  Read raw data file header from file
#  fd:  raw data file descriptor
def Read_FHeader(df):
  p = PF(df)		#read header data

  fs = float(p.get("adc_samp_freq")[0])
  ns = int(p.get("adc_num_samp")[0])
  psf = int((p.get("adc_presum")[0]))
  ncycle = int(p.get("num_cycles")[0])
  cycle_delay = float(p.get("cycle_delay")[0])
  elev0 = float(p.get("start_elev_ang")[0])
  elev_step = float(p.get("elev_step")[0])
  nelev = int(p.get("elev_num_step")[0])
  pmodes = p.get("measure_modes")	# string of modes HL,CL,AL,HP,VP

  hdr_len = struct.calcsize('=iiiffffffffi') #length of data record header
  rec_len = ns/psf * 8           #number of bytes in a data record
  return (fs, ns, psf, ncycle, cycle_delay, elev0, elev_step, nelev, pmodes, hdr_len, rec_len)

## (year,month,day,sod,elev_ang,t_plate,t_sink,t_ext,t_cable,dc_offset1,dc_offset2,nps) = Read_LHeader(df)
#  Read raw data line header from file
#
#  Read raw data line header from file
#  fd:  raw data file descriptor
def Read_LHeader(df):
  hdr_len = struct.calcsize('=iiiffffffffi') #length of data record header
  hdr = df.read(hdr_len)
  if len(hdr) < hdr_len:
    raise EOF, 'EOF'

  hdr_data= struct.unpack('=iiiffffffffi',hdr)
  year = hdr_data[0] 
  month = hdr_data[1] 
  day  = hdr_data[2]
  sod  = hdr_data[3]
  elev_ang = hdr_data[4]
  t_plate = hdr_data[5]
  t_sink = hdr_data[6]
  t_ext = hdr_data[7]
  t_cable = hdr_data[8]
  dc_offset1 = hdr_data[9]
  dc_offset2 = hdr_data[10]
  nps  = hdr_data[11]
  return (year,month,day,sod,elev_ang,t_plate,t_sink,t_ext,t_cable,dc_offset1,dc_offset2,nps)


## (h,m,s) = Sod2hms(sod)
#  Get hours,minutes,seconds from seconds of day
#
#  Get hours,minutes,seconds from seconds of day
#  sod:  seconds of day
def Sod2hms(sod):
  h,s = divmod(sod,3600)
  m,s = divmod(s,60)
  return (h,m,s)


