"""Check a proposed public ZIP without opening any credentialed database.
Usage: python tests/check_release.py package.zip [package2.zip ...]
Only standard-library modules are used. Error messages never print data values.
"""
from pathlib import Path, PurePosixPath
import csv, io, re, sys, zipfile, json

ID = re.compile(r'(?i)(^|_)(subject_?id|hadm_?id|stay_?id|patientunitstayid|uniquepid|hospital_?id|healthsystemstayid|patient_?id|admittime|dischtime|deathtime|charttime|panel_time|dateofbirth|dob)($|_)')
PRIVATE_NAMES = {'db.env', '.pgpass', 'pgpass.conf', '.env', '.rdata', '.rhistory'}
ALLOWLIST = Path(__file__).with_name('release_csv_allowlist.json')
def audit(path):
    bad=[]; csvs=0
    policy = json.loads(ALLOWLIST.read_text(encoding='utf-8'))
    declared = policy.get(Path(path).name, {})
    if not declared:
        m = re.search(r"(?i)additional[_ ]?file[_ ]?(\d)", Path(path).name)
        if m:
            for k, v in policy.items():
                if re.search(r"(?i)additional[_ ]?file[_ ]?" + m.group(1) + r"($|[^0-9])", k):
                    declared = v
                    break
    with zipfile.ZipFile(path) as z:
        if z.testzip() is not None:bad.append('ZIP CRC failure')
        for info in z.infolist():
            if info.is_dir():continue
            name=info.filename.replace('\\','/');p=PurePosixPath(name)
            if '..' in p.parts or name.startswith('/') or ':' in name:bad.append('unsafe archive path: '+name)
            if p.name.lower() in PRIVATE_NAMES or any(part in {'interim','logs','Rlib','.git'} for part in p.parts):bad.append('restricted file/path: '+name)
            if p.name.startswith('roc_coords') or p.name in {'eicu_sicu_aligned_panels.csv','eicu_sicu_patient_index_roster.csv','eicu_sicu_external_validation_dataset.csv'}:bad.append('restricted derived data: '+name)
            if p.suffix.lower() in {'.rds','.rdata','.parquet','.feather','.sqlite','.db','.gz','.zip'}:bad.append('binary/nested data require separate review: '+name)
            if p.suffix.lower()=='.csv':
                csvs+=1
                allowed = declared.get(name)
                if allowed is None:
                    allowed = declared.get(name.split("/", 1)[-1])
                if allowed is None:bad.append('undeclared CSV: '+name)
                reader=csv.DictReader(io.StringIO(z.read(info).decode('utf-8-sig')))
                fields=reader.fieldnames or []
                if not fields or any(ID.search(f) for f in fields):bad.append('identifier/time column in CSV: '+name)
                rows = list(reader)
                if allowed and (fields != allowed['columns'] or len(rows) > allowed['max_rows']):bad.append('CSV exceeds reviewed schema/row bound: '+name)
                if len(fields) != len(set(fields)):bad.append('duplicate CSV columns: '+name)
                if any(None in row or any(v is None for v in row.values()) for row in rows):bad.append('nonrectangular CSV: '+name)
    return bad,csvs
if __name__=='__main__':
    if len(sys.argv)<2:raise SystemExit(__doc__)
    if not ALLOWLIST.exists():raise SystemExit('release CSV allowlist is missing')
    failures=0
    for arg in sys.argv[1:]:
        issues,n=audit(arg);failures+=len(issues)
        print(f'{Path(arg).name}: {"FAIL" if issues else "PASS"}; {n} CSVs')
        for issue in issues:print('  '+issue)
    raise SystemExit(1 if failures else 0)
