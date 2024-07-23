# functions related to hdf5 files
import h5py
import numpy as np

def read_IC_hdf5(filename):
 data = {}
 file = h5py.File(filename,'r')
 keys = list(file.keys())
 for key in keys:
  if key == 'Header':
   header = dict(file['Header'].attrs.items())
   data['Header'] = header
  else:
   partData = file.get(key)
   partKeys = list(partData.keys())
   data[key] = {}
   for partKey in partKeys:
    data[key][partKey] = np.array(partData.get(partKey))
 file.close()
 return data

def write_IC_hdf5(filename, data):
 def ifset(d,key):
  if key in d.keys():
   result=d[key]
  else:
   result=None
   if key in ['lt','fmt']:  #for plot
    result=''
  return result

 def ifset2(d,key,value):
  if value is None:
   result=ifset(d,key)
  else:
   result=value
  return result

 if isinstance(data, dict):
  data=[data]

 BoxSize = None
 NumPartType = 6
 MassTable = np.zeros(NumPartType, dtype = float)
 NumPart = np.zeros(NumPartType, dtype = int)
 
 for d in data:
  BoxSize = ifset2(d, 'BoxSize', BoxSize)
  i = d['PartType']
  MassTable[i] = d['PartMass']
  NumPart[i] = d['count']
  
 file = h5py.File(filename+'.hdf5','w')

 group = file.create_group("Header")
 if BoxSize is not None:
  group.attrs.create("BoxSize", BoxSize, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 else:
  group.attrs.create("BoxSize", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Flag_Cooling", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_DoublePrecision", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_Feedback", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_IC_Info", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_Metals", 0, shape=None, dtype=h5py.h5t.STD_I32LE) #in makegal ics is 1
 group.attrs.create("Flag_Sfr", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("Flag_StellarAge", 0, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("HubbleParam", 1, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("MassTable", MassTable, shape=None,
                                                     dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("NumFilesPerSnapshot", 1, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("NumPart_ThisFile", NumPart, shape=None, dtype=h5py.h5t.STD_I32LE)
 group.attrs.create("NumPart_Total", NumPart, shape=None, dtype=h5py.h5t.STD_U32LE)
 group.attrs.create("NumPart_Total_HighWord", (0,0,0,0,0,0), shape=None, dtype=h5py.h5t.STD_U32LE)
 group.attrs.create("Omega0", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("OmegaLambda", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Redshift", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)
 group.attrs.create("Time", 0, shape=None, dtype=h5py.h5t.IEEE_F64LE)

 ID_offset = 0
 for d in data:
  group = file.create_group("PartType"+str(d['PartType']))
  dataset = group.create_dataset("Coordinates", (d['count'],3), data=d['Coordinates'],
                                                            dtype=h5py.h5t.IEEE_F32LE)
  dataset = group.create_dataset("ParticleIDs", (d['count'],),
                   data=np.array(range(ID_offset,ID_offset+d['count'])), dtype=h5py.h5t.STD_I32BE)
  ID_offset += d['count']
  dataset = group.create_dataset("Velocities", (d['count'],3), data=d['Velocities'],
                                                            dtype=h5py.h5t.IEEE_F32LE)
  if d["PartType"] == 0:
   dataset = group.create_dataset("InternalEnergy", (d['count'],),
                   data=d["InternalEnergy"], dtype=h5py.h5t.IEEE_F32LE)
 file.close()