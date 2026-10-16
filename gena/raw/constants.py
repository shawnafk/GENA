import os

from common.constants import MASS_E, MASS_H, MASS_O, SPEED_C
from gena.constants import init_param

#----------------------------------------------- for decoding 
#packet num
#os.environ["GENA_CHUNK_SIZE"] = str(int(86400*32*20))
if 'GENA_CHUNK_SIZE' not in os.environ:
    chunk_size = int(100000)
else:
    chunk_size = int(os.environ['GENA_CHUNK_SIZE'])

if 'GENA_USE_REPAIR' not in os.environ:
    repair = 0
else:
    repair = int(os.environ['GENA_USE_REPAIR'])

if 'GENA_SAVE_L1B' not in os.environ:
    #default is save l1b
    save_l1b = 1
else:
    save_l1b = int(os.environ['GENA_SAVE_L1B'])

if 'GENA_REPAIR_METHOD' not in os.environ:
    method = 2
else:
    method = int(os.environ['GENA_REPAIR_METHOD'])

if 'GENA_PDF_ONLY' not in os.environ:
    #default is convert all
    pdf_only = 0
else:
    pdf_only = int(os.environ['GENA_PDF_ONLY'])

if 'GENA_SAVE_PDF' not in os.environ:
    #default is save pdf
    save_pdf = 1
else:
    save_pdf = int(os.environ['GENA_SAVE_PDF'])

if 'GENA_SAVE_L1' not in os.environ:
    #default is not save l1
    save_l1 = 0
else:
    save_l1 = int(os.environ['GENA_SAVE_L1'])

if 'GENA_SAVE_L1A' not in os.environ:
    #default is save l1a
    save_l1a = 1
else:
    save_l1a = int(os.environ['GENA_SAVE_L1A'])

if 'GENA_SAVE_L2' not in os.environ:
    #default is save l2
    save_l2 = 1
else:
    save_l2 = int(os.environ['GENA_SAVE_L2'])



if __name__ == '__main__':
    print('test')
