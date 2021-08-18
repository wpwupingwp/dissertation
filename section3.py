#!/usr/bin/python3

from multiprocessing import Pool
from pathlib import Path
from subprocess import run
from sys import argv
import logging


FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
TEMP_LOG = 'align_to_ref.log'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO,
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler(TEMP_LOG)])
import coloredlogs
coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
log = logging.getLogger()


def copy_():
    folder = Path('./raw-newest')
    dest = Path(Path(argv[1]).stem)
    bad = []
    with open(argv[1], 'r') as _:
        for line in _:
            try:
                f = list(folder.glob(line.strip()+'*'))[0]
                cmd = f'cat "{f}" >> {dest/dest.name}.gb'
                print(cmd)
                r = run(cmd, shell=True)
                if r.returncode != 0:
                    bad.append(cmd)

            except IndexError:
                pass
    return bad


def run_cmd(f):
    folder = Path(f)
    gb = folder / Path(f+'.gb')
    cmd = (f'python3 -m BarcodeFinder.gb2fasta -gb {gb} '
           f'-out {folder}_divide -rename -unique no '
           f'-allow_repeat -allow_invert_repeat')
    print(cmd)
    r = run(cmd, shell=True)
    print(r)
    return r.returncode


def divide():
    folders = ('all_species', 'all_family', 'all_genus', 'all_sample')
    pool = Pool()
    with Pool() as pool:
        pool.map(run_cmd, folders)
        # apply_async()
        # pool.close()
        # pool.join()
    print('divide done')


def align():
    folders = ('all_family', 'all_genus', 'all_species', 'all_sample')
    new_folders = [Path(i+'_divide')/'Unique' for i in folders]
    skip = ('gene-rps12', )
    for folder in new_folders:
        log.info(f'Folder {folder}')
        fasta = list(folder.glob('*.uniq'))
        log.info(f'{len(fasta)} files')
        bad = open(folder.stem+'_bad.txt', 'w')
        n = 0
        n_bad = 0
        for fas in fasta:
            n += 1
            log.info(fas.name)
            if skip[0] in fas.name:
                log.warning(f'Skip {fas}')
            cmd = (f'mafft --reorder --adjustdirection --thread 15 '
                   f'--genafpair --maxiterate 1000 '
                   f'{fas} > {fas.with_suffix(".aln")}')
            r = run(cmd, shell=True)
            if r.returncode != 0:
                bad.write(str(fas)+'\n')
                log.warning(f'{fas} failed.')
                n_bad += 1
        log.info(f'Total {n}')
        log.info(f'Bad {n_bad}')


def align_to_ref():
    folders = ('all_genus', 'all_species', 'all_sample')
    ref_ = Path('all_family_divide') / 'Alignment'
    ref = list(ref_.glob('*.fasta'))
    ref_dict = {i.stem: i for i in ref}

    new_folders = [Path(i+'_divide')/'Unique' for i in folders]
    skip = ('gene-rps12', 'CDS-ycf1', 'CDS-ycf2')
    # hard to align
    # hard = ('gene-ycf1', 'gene-ycf2')
    for folder in new_folders:
        log.info(f'Folder {folder}')
        if 'all_genus' in str(folder.name):
            fasta = list(folder.glob('*.uniq'))
        else:
            fasta = list(folder.glob('spacer*.uniq'))
        log.info(f'{len(fasta)} files')
        bad = open(folder.stem+'_bad.txt', 'w')
        n = 0
        n_bad = 0
        for fas in fasta:
            out = fas.with_suffix('.aln')
            if out.exists():
                log.warning(f'Skip existing file {out}')
                continue
            n += 1
            log.info(fas.name)
            if fas.stem in skip:
                log.warning(f'Skip {fas}')
                continue
            if fas.stem in ref_dict:
                ref_fasta = ref_dict[fas.stem]
                log.info(f'Use ref {ref_fasta} to align {fas}')
                cmd = (f'mafft --add {fas} --reorder --adjustdirection '
                       # f'--thread 15 --genafpair --maxiterate 1000 '
                       f'--auto --thread 15'
                       f'{ref_fasta} > {fas.with_suffix(".aln")}')
            else:
                cmd = (f'mafft --reorder --adjustdirection --thread 15 '
                       f'--auto --thread 15'
                       f'{fas} > {fas.with_suffix(".aln")}')
            r = run(cmd, shell=True)
            if r.returncode != 0:
                bad.write(str(fas)+'\n')
                log.warning(f'{fas} failed.')
                n_bad += 1
        log.info(f'Total {n}')
        log.info(f'Bad {n_bad}')


def run_evaluate(f):
    cmd = (f'python3 -m BarcodeFinder.evaluate -quick -aln {f} '
           f'-out ~/analyze/all_genus_divide/Evaluate/{f.stem}-out')
    print(cmd)
    r = run(cmd, shell=True)
    return r.returncode


def evaluate():
    folder = Path('~/analyze/all_genus_divide/Alignment').expanduser()
    ref = Path('~/analyze/all_family_divide/Alignment').expanduser()
    ref_fasta = ref.glob('*.fasta')
    ref_set = set([i.name for i in ref_fasta])
    # skip ycf1 and ycf2
    files = [i for i in folder.glob('*.fasta') if i.name in ref_set]
    with Pool(15) as pool:
        pool.map(run_evaluate, files)

evaluate()
