'''Create the anonymous submission archive, excluding survey microdata.'''
from pathlib import Path
import hashlib
import json
import os
import re
import sys
from zipfile import ZipFile,ZipInfo,ZIP_DEFLATED

sys.dont_write_bytecode=True
P=Path(__file__).resolve().parent
ARCHIVE=P.parent/'Anonymous_peer_review_code.zip'
EXCLUDE={'checksums.sha256','package_audit.json',ARCHIVE.name}

def included(f):
    relative=f.relative_to(P)
    return (f.is_file() and relative.parts[0]!='local_data' and
        '__pycache__' not in relative.parts and f.suffix not in ['.pyc','.zip'] and relative.name not in EXCLUDE)

def main():
    files=sorted(f for f in P.rglob('*') if included(f))
    findings=[]
    personal=os.environ.get('USERNAME','')
    pattern=re.compile(r'(?i)[A-Z]:[\\/](?:Users|Second_Class)[\\/]')
    for f in files:
        if f.suffix.lower() in {'.py','.r','.ps1','.cmd','.txt','.csv','.json','.md','.gmt'}:
            text=f.read_text(encoding='utf-8-sig')
            if pattern.search(text):findings.append({'file':f.relative_to(P).as_posix(),'type':'personal_absolute_path'})
            if len(personal)>3 and personal.casefold() in text.casefold():
                findings.append({'file':f.relative_to(P).as_posix(),'type':'local_account_name'})
        if f.suffix=='.pdf':
            import pymupdf
            with pymupdf.open(f) as doc:
                for key in ['author','creator','subject','title']:
                    value=doc.metadata.get(key,'') or ''
                    if pattern.search(value) or (len(personal)>3 and personal.casefold() in value.casefold()):
                        findings.append({'file':f.relative_to(P).as_posix(),'type':'pdf_identity_metadata'})
    if findings:raise RuntimeError('Anonymity audit found: '+json.dumps(findings))
    # Include names and hashes only, never absolute source paths.
    checksums=''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.relative_to(P).as_posix()+'\n' for f in files)
    (P/'checksums.sha256').write_text(checksums,encoding='utf-8')
    files.append(P/'checksums.sha256')
    count=sum(f.is_file() for f in P.rglob('*'))
    count+=int(not (P/'package_audit.json').exists())
    if count>100:raise RuntimeError(f'Recursive file count {count} exceeds 100.')
    directories={'.':count}
    for folder in sorted(f for f in P.rglob('*') if f.is_dir()):
        directories[folder.relative_to(P).as_posix()]=sum(f.is_file() for f in folder.rglob('*'))
    report={'status':'passed','physical_files_in_code':count,'limit_per_directory_recursive':100,
        'directory_recursive_counts':directories,'submission_files':len(files)+1,
        'anonymity_text_and_pdf_scan':'passed','local_survey_microdata_in_submission':False,
        'archive_member_timestamps':'fixed','symlinks_allowed':False,
        'binary_inputs':'plain expression/design lists and count arrays; no fitted model environments',
        'data_access':'See DATA_ACCESS.md; complete local data are excluded from submission.'}
    (P/'package_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    files.append(P/'package_audit.json')
    if any(f.is_symlink() for f in files):raise RuntimeError('Symlinks are not allowed in the submission.')
    with ZipFile(ARCHIVE,'w',compression=ZIP_DEFLATED,compresslevel=6) as z:
        for f in files:
            info=ZipInfo(f.relative_to(P).as_posix(),date_time=(1980,1,1,0,0,0))
            info.compress_type=ZIP_DEFLATED;info.external_attr=0o644<<16
            z.writestr(info,f.read_bytes())
    with ZipFile(ARCHIVE) as z:
        if z.testzip():raise RuntimeError('Archive CRC verification failed.')
        if any(n.startswith('local_data/') for n in z.namelist()):raise RuntimeError('Local data entered submission.')
        if len(z.namelist())>100:raise RuntimeError('Submission archive contains more than 100 files.')
    print(json.dumps(report,indent=2));print('Created in parent directory:',ARCHIVE.name)

if __name__=='__main__':main()
