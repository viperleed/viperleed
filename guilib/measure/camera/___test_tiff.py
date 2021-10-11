
import os, sys
import timeit
from time import perf_counter as timer
import zipfile
import subprocess

import numpy as np
import cv2
import h5py

from scipy import stats

base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tmp_files")

infile = "in2.tif"

# NB: all field values are left-justified in the 4 bytes at the end of the tags

tag_ids = {254: 'NewSubfileType',             # NO
           255: 'SubfileType',                # NO
           256: 'ImageWidth',                 # b'\x01\x00' + b'\x00\x04' (may also be b'\x00\x03') + b'\x00\x00\x00\x01' + bytes(n_cols)
           257: 'ImageLength',  # = height    # b'\x01\x01' + b'\x00\x04' (may also be b'\x00\x03') + b'\x00\x00\x00\x01' + bytes(n_rows)
           258: 'BitsPerSample',              # b'\x01\x02\x00\x03\x00\x00\x00\x01' + bytes(n_bits)
           259: 'Compression',                # NO
           262: 'PhotometricInterpretation',  # b'\x01\x06\x00\x03\x00\x00\x00\x01\x00\x01\x00\x00' (black is zero)
           263: 'Thresholding',               # NO
           269: 'DocumentName',               # NO
           270: 'ImageDescription',           # NO, perhaps later
           271: 'Make',                       # NO, perhaps later
           272: 'Model',                      # NO, perhaps later
           273: 'StripOffsets',               # YES: where do the data strips start (wrt beginning of file)
                                              #      b'\x01\x11' + (word [b'\x00\x03'] or dword [b'\x00\x04']) + (depends on n_strips, should be b'\x00\x00\x00\x01' if we don't split images) + bytes(offset address)
                                              #      offset address needs to be contiguous when no compression, at the end of the tags
           274: 'Orientation',                # NO (defaults to 1st row = top, 1st col = left)
           277: 'SamplePerPixel',             # NO
           278: 'RowsPerStrip',               # in principle not necessary, although it should be < 8K if one wants to speed up access to the file. e.g., n_rows can be used:
                                              # b'\x01\x16' + b'\x00\x03' (may also be b'\x00\x04') + b'\x00\x00\x00\x01' + bytes(n_rows)
           279: 'StripByteCounts',            # b'\x01\x17' + b'\x00\x04' (may also be b'\x00\x03') + b'\x00\x00\x00\x01' + bytes(n_cols*rows_per_strip(=n_rows)*n_bits/8)
           282: 'XResolution',                # NO
           283: 'YResolution',                # NO
           284: 'PlanarConfiguration',        # NO
           286: 'XPosition',                  # NO
           287: 'YPosition',                  # NO
           290: 'GrayResponseUnit',           # NO
           291: 'GrayResponseCurve',          # NO
           292: 'Group3Options',              # NO
           293: 'Group4Options',              # NO
           296: 'ResolutionUnit',             # NO
           297: 'PageNumber',                 # NO, perhaps later for energy?
           301: 'ColorResponseCurves',        # NO
           305: 'Software',                   # NO, perhaps later
           306: 'DateTime',                   # NO, perhaps later
           315: 'Artist',                     # NO
           316: 'HostComputer',               # NO, perhaps later
           317: 'Predictor',                  # NO
           318: 'WhitePoint/ColorImageType',  # NO - depends on whether it is WORD(3)/RATIONAL(5)
           319: 'PrimaryChromaticities/ColorList', # NO - depends if it is RATIONAL(5)/BYTE(1)orWORD(3)
           320: 'ColorMap',                   # NO
           }

tag_bytes = {1: 1,
             2: 1, # 'ascii', last byte is b'\x00'
             3: 2,
             4: 4,
             5: 8
             }

def read_uncompressed_tiff(fname):
    # TIFF header:
    # [0:2] 2 chars (w = 2), 'II' (-> 'little') or 'MM' (-> 'big')
    # [2:4] 1 word  (w = 2), tiff version, typically = 42
    # [4:8] 1 dword (w = 4), offset of image start
    #
    # From image start offset
    # [0:2] number of entries, each entry is 12 bytes long
    # 
    # Each entry:
    # [0:2] id
    # [2:4] type
    # [4:8] length
    # [8:12] value (if <= 4bytes) or address to value
    #
    # write a dword zero after the end of the tags (at tagstart + 2 + numtags*12)
    
    
    f = open(fname, 'rb')
    f.seek(0)
    b = f.read(8)
    byteorder = {b'II': 'little', b'MM': 'big'}[b[:2]]
    version = int.from_bytes(b[2:4], byteorder)
    tag_start = int.from_bytes(b[4:], byteorder)
    
    f.seek(tag_start)
    b = f.read(2)
    num_entries = int.from_bytes(b[:2], byteorder)
    
    tags = {}
    # read tags, each is 12 bytes long
    for i in range(num_entries):
        b = f.read(12)
        # print(b)
        tag = {}
        id = int.from_bytes(b[:2], byteorder)
        typ = int.from_bytes(b[2:4], byteorder)
        tag['id'] = tag_ids.get(id, id)
        tag['size'] = tag_bytes[typ]
        tag['length'] = int.from_bytes(b[4:8], byteorder)
        tag['offset'] = int.from_bytes(b[8:], byteorder)
        
        size = tag['length']*tag['size']
        if size <= 4:
            value_bytes = b[8:8+size]
        else:
            start = tag['offset']
            curpos = f.tell()
            f.seek(start)
            str = ''
            value_bytes = f.read(size)
            f.seek(curpos)
        # print(value_bytes)
        if typ == 2:  # string
            value = value_bytes[:-1].decode('utf-8')
        else:  # read into a list
            vals = (value_bytes[i:i+tag['size']]
                    for i in range(0, size, tag['size']))
            if typ == 5:  # ratio
                value = tuple(int.from_bytes(v[:4], byteorder)
                              / int.from_bytes(v[4:], byteorder)
                              for v in vals)
            else:
                value = tuple(int.from_bytes(v, byteorder) for v in vals)
            if (len(value) == 1
                    and not tag['id'] in ('StripOffsets', 'StripByteCounts')):
                value = value[0]
        tag['value'] = value
        tags.update({tag['id']: tag['value']})
    
    # print(tags)
    
    # finally, get the data: we expect ImageWidth * ImageLength entries of
    # BitsPerSample/8 length (in bytes), arranged into strips of RowsPerStrip
    # number of rows, each starting at one of the StripOffsets
    data = []
    bytes_per_point = int(tags['BitsPerSample']/8)
    for b_count, offset in zip(tags['StripByteCounts'], tags['StripOffsets']):
        f.seek(offset)
        b = f.read(b_count)  # todo: last strip may have less data
        data.extend(int.from_bytes(b[i:i+bytes_per_point], byteorder)
                    for i in range(0, b_count, bytes_per_point))
    f.close()
    
    data = np.asarray(data, dtype='>u2').reshape(tags['ImageLength'],
                                                 tags['ImageWidth'])
    if tags['PhotometricInterpretation'] == 0:
        data = 2**tags['BitsPerSample'] - 1 - data
    
    # print(data)
    return data

tiff_tags = (
    b'MM\x00\x2a\x00\x00\x00\x08\x00\x07'                  # header
    + b'\x01\x00\x00\x04\x00\x00\x00\x01%(n_cols)b'        # 'ImageWidth'
    + b'\x01\x01\x00\x04\x00\x00\x00\x01%(n_rows)b'        # 'ImageLength'
    + b'\x01\x02\x00\x03\x00\x00\x00\x01\x00\x10\x00\x00'  # 'BitsPerSample'
    + b'\x01\x06\x00\x03\x00\x00\x00\x01\x00\x01\x00\x00'  # no.262
    + b'\x01\x11\x00\x04\x00\x00\x00\x01\x00\x00\x00\x62'  # 'StripOffsets'
    + b'\x01\x16\x00\x04\x00\x00\x00\x01%(n_rows)b'        # 'RowsPerStrip'
    + b'\x01\x17\x00\x04\x00\x00\x00\x01%(b_count)b'       # 'StripByteCounts'
    + b'\x00\x00\x00\x00'                                  # IFD termination
    )


def write_basic_tiff(fname, data):
    """
    Write the np.array data to fname as a 16-bit TIFF file
    """
    data = np.asarray(data, dtype='>u2')
    n_rows, n_cols = data.shape
    
    # header:  'MM'(2) '42'(2) start_of_tags(4) 'no.tags[=7]'(2)
    # header = b'MM\x00\x2a\x00\x00\x00\x08\x00\x07'
    # tags = b''
    
    # The seven tags are:
    # (1) 'ImageWidth' (no.256)
    # tags += b'\x01\x00\x00\x04\x00\x00\x00\x01' + n_cols.to_bytes(4, 'big')
    # (2) 'ImageLength', i.e. height (no.257)
    # tags += b'\x01\x01\x00\x04\x00\x00\x00\x01' + n_rows.to_bytes(4, 'big')
    # (3) 'BitsPerSample' (no.258); 16-bit, left-justified to the first word
    # tags += b'\x01\x02\x00\x03\x00\x00\x00\x01' + b'\x00\x10\x00\x00'
    # (4) 'PhotometricInterpretation' (no.262); black is zero
    # tags += b'\x01\x06\x00\x03\x00\x00\x00\x01\x00\x01\x00\x00'
    # (5) 'StripOffsets' (n.273); will have one strip only
    # datastart inlcudes also 4 = len(b'\x00\x00\x00\x00'), due to the end of
    # the IFD. Thus
    #     datastart = len(header) + 12*n_tags + 4 = 10 + 12*7 + 4 = 98
    # i.e., datastart = b'\x00\x00\x00\x62'
    # tags += b'\x01\x11\x00\x04\x00\x00\x00\x01\x00\x00\x00\x62'
    # (6) 'RowsPerStrip' (no.278)
    # tags += b'\x01\x16\x00\x04\x00\x00\x00\x01' + n_rows.to_bytes(4, 'big')
    # (7) 'StripByteCounts' (no.279)
    # tags += (b'\x01\x17\x00\x04\x00\x00\x00\x01'
             # + int(n_cols * n_rows * 2).to_bytes(4, 'big'))
    
    tags = tiff_tags %{b'n_cols': n_cols.to_bytes(4, 'big'),
                       b'n_rows': n_rows.to_bytes(4, 'big'),
                       b'b_count': (n_cols * n_rows * 2).to_bytes(4, 'big')}
    try:
        with open(fname, 'wb') as f:
            f.write(tags + data.tobytes())
    except:
        pass


def write_hd5(fname, data):
    h5_kwargs = {'chunks': True,
                 'compression': 'gzip',
                 'compression_opts': 4}
    
    data = np.asarray(data, dtype='>u2')
    
    f = h5py.File(fname, 'w')
    dset = f.create_dataset("img", data=data, **h5_kwargs)
    f.close()

def read_multi_tiff(dir_name, maxread=None):
    files = [f for f in os.listdir(dir_name) if ".tif" in f]
    
    if maxread is None:
        maxread = len(files)
    files = files[0:maxread]
    
    # return np.asarray([cv2.imread(os.path.join(dir_name, f),
                                  # cv2.IMREAD_GRAYSCALE) for f in files])
    
    return np.asarray([read_uncompressed_tiff(os.path.join(dir_name, f))
                       for f in files], dtype='>u2')
    
# read_uncompressed_tiff(os.path.join(base_path, 'in.tif'))
# d = read_uncompressed_tiff(os.path.join(base_path, infile))

outfile = 'out.tif'
fname = os.path.join(base_path, outfile)
fname_cv2 = os.path.join(base_path, 'cv2out.tif')
fname_h5 = os.path.join(base_path, 'h5out.h5')
tmp_dir = os.path.join(base_path, 'tmp')
tmp_dir_in = os.path.join(base_path, 'tmp_in')
tmp_dir_out = os.path.join(base_path, 'tmp_out')
archive_name = os.path.join(base_path, 'out.zip')
archive_cv_name = os.path.join(base_path, 'out_cv.zip')
archive_zip_name = os.path.join(base_path, 'out_zip.zip')
path_7zip = r"E:\Software\Program Files\7-Zip\7z.exe"

def multi_tiffs(archive_name, tmp_dir, data, repeat=900):
    for i in range(repeat):
        write_basic_tiff(os.path.join(tmp_dir, f'out{i}.tif'), data)
    pack_zip(archive_name, tmp_dir, contains='out')

def multi_cv2_tiffs(archive_name, tmp_dir, data, repeat=900):
    for i in range(repeat):
        cv2.imwrite(os.path.join(tmp_dir, f'cv2{i}.tif'), data)
    pack_zip(archive_name, tmp_dir, contains='out')

def multi_h5(fname, data, repeat=900):
    h5_kwargs = {'chunks': True,
                 'compression': 'gzip',
                 'compression_opts': 4}
    
    data = np.asarray([np.asarray(data, dtype='>u2') for i in range(repeat)])
    
    with h5py.File(fname, 'w') as f:
        f.create_dataset("img", data=data, **h5_kwargs)

def multi_h5_dlist(fname, data):
    h5_kwargs = {'chunks': True,
                 'compression': 'gzip',
                 'compression_opts': 4}
    
    with h5py.File(fname, 'w') as f:
        f.create_dataset("img", data=data, **h5_kwargs)

def pack_zip(archive_name, files_dir, contains='', use_7z=True):
    if not use_7z:
        with zipfile.ZipFile(archive_name,
                             'w',
                             compression=zipfile.ZIP_DEFLATED,
                             compresslevel=4) as archive:
            files = [f
                     for f in os.listdir(files_dir)
                     if (".tif" in f) and (contains in f)]
            for f in files:
                f_path = os.path.join(files_dir, f)
                # add files to archive
                archive.write(f_path, os.path.basename(f_path))
                
                # # and remove them
                # try:
                    # os.unlink(f_path)
                # except Exception as exc:
                    # print('Failed to delete %s. Reason: %s' % (f_path, exc))
        return
    os.chdir(files_dir)
    zip_proc = subprocess.check_output([path_7zip, "a", "-tzip", archive_name, f"{contains}*"])
    os.chdir(base_path)

def append_to_zip(archive_name, files_dir, contains='',
                  use_7z=True, maxfiles=None):
    files = [f
             for f in os.listdir(files_dir)
             if all(c in f for c in (".tif", contains))]
    if maxfiles is None:
        maxfiles = len(files)
    files = files[:maxfiles]
    
    dt = []
    for f in files:
        t0 = timer()
        if not use_7z:
            with zipfile.ZipFile(archive_name,
                                 'a',
                                 compression=zipfile.ZIP_DEFLATED,
                                 compresslevel=9) as archive:
                f_path = os.path.join(files_dir, f)
                archive.write(f_path, os.path.basename(f_path))
        else:
            os.chdir(files_dir)
            zip_proc = subprocess.check_output([path_7zip, "a", "-tzip", archive_name, str(f)])
            os.chdir(base_path)
        t1 = timer()
        dt.append(t1-t0)
        print(t1-t0, end=' ')
    print()
    return np.asarray(dt)

def time_me(nruns=1000):
    # t_mine = timeit.Timer(stmt="tst.write_basic_tiff(tst.fname, tst.d)",
                          # setup="""
# import test_tiff as tst
# tiff_tags = tst.tiff_tags""")
    
    # t_cv2 = timeit.Timer(stmt="tst.cv2.imwrite(tst.fname_cv2, tst.d)",
                         # setup="import test_tiff as tst")
    
    # t_h5 = timeit.Timer(stmt="tst.write_hd5(tst.fname_h5, tst.d)",
                         # setup="import test_tiff as tst")
    
    # print(f"Mine: {t_mine.timeit(nruns)*1000/nruns} ms")
    # print(f"cv2: {t_cv2.timeit(nruns)*1000/nruns} ms")
    # print(f"hd5: {t_h5.timeit(nruns)*1000/nruns} ms")
    
    print("----------------- MULTI ------------------")
    
    # t_multi_mine = timeit.Timer(
        # stmt="tst.multi_tiffs(tst.archive_name, tst.tmp_dir, tst.d)",
        # setup="""
# import test_tiff as tst
# tiff_tags = tst.tiff_tags""")
    

    # t_multi_cv = timeit.Timer(
        # stmt="tst.multi_cv2_tiffs(tst.archive_cv_name, tst.tmp_dir, tst.d)",
        # setup="import test_tiff as tst")
    
    t_multi_h5 = timeit.Timer(
        # stmt="tst.multi_h5(tst.fname_h5, tst.d)",
        stmt="tst.multi_h5_dlist(tst.fname_h5, tst.d)",
        setup="import test_tiff as tst")
    
    t_pack_zip = timeit.Timer(
        stmt="tst.pack_zip(tst.archive_zip_name, tst.tmp_dir_out)",
        setup="import test_tiff as tst")
    
    # print(f"Mine: {t_multi_mine.timeit(nruns)*1000/nruns} ms")
    # print(f"cv2: {t_multi_cv.timeit(nruns)*1000/nruns} ms")
    print(f"hd5: {t_multi_h5.timeit(nruns)*1000/nruns} ms")
    print(f"7z: {t_pack_zip.timeit(nruns)*1000/nruns} ms")
    

if __name__ == '__main__':
    # write_basic_tiff(os.path.join(base_path, outfile), d)
    # cv2.imwrite(os.path.join(base_path, 'cv2out.tif'), d)
    
    # time_me(1)
    # multi_tiffs(archive_name, tmp_dir, d, repeat=10)
    # multi_cv2_tiffs(archive_cv_name, tmp_dir, d, repeat=10)
    # multi_h5(fname_h5, d)
    
    # d = read_multi_tiff(tmp_dir_in)
    # for i, di in enumerate(d):
        # write_basic_tiff(os.path.join(tmp_dir_out, f'{i}.tif'), di)
    
    # print('\n RE-READ\n')
    # read_uncompressed_tiff(os.path.join(base_path, outfile))
    
    # t0=timer()
    # multi_h5_dlist(fname_h5, d)
    # print(f"hd5: {timer()-t0} s")
    
    # t0=timer()
    # pack_zip(archive_zip_name, tmp_dir_in)
    # print(f"7z: {timer()-t0} s")
    MF = 881
    
    dt_zip = append_to_zip(archive_zip_name, tmp_dir_in, maxfiles=MF)
    dt_py = append_to_zip(archive_name, tmp_dir_in, use_7z=False, maxfiles=MF)
    
    print()
    print("7z:", stats.describe(dt_zip))
    print("zipfile:", stats.describe(dt_py))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    