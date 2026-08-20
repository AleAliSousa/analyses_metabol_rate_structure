"""Python port of compare_stephan.R's logic, run against the local snapshots.

R is unavailable in the sandbox and OneDrive is not mounted, so this exercises
the same alias / pairing / comparison / source-join logic on real repo data
plus a merge-form mock of volumes_wide_select.csv, to catch logic errors before
Alexandra runs the R version against live OneDrive.
"""
import csv, re, sys, random
from collections import Counter, defaultdict

ROOT = '/sessions/nifty-gracious-darwin/mnt/analyses_metabol_rate_structure/'
FAIL = []
def check(name, cond, detail=''):
    print(('  PASS  ' if cond else '  FAIL  ') + name + (('  -- ' + str(detail)) if detail else ''))
    if not cond: FAIL.append(name)

# ---------- port of R/species_aliases.R --------------------------------------
ALIAS = [  # tree_label, merge_label, kind
    ("Callithrix_pygmaea",   "Callithrix pygmaea",    "spelling"),
    ("Lagothrix_lagotricha", "Lagothrix lagothricha", "spelling"),
    ("Plecturocebus_moloch", "Callicebus moloch",     "genus"),
    ("Gorilla_gorilla",      "Gorilla sp.",           "pooled"),
    ("Tarsius_syrichta",     "Tarsius sp.",           "pooled"),
    ("Pongo_abelii",         "Pongo sp.",             "ambiguous"),
    ("Pongo_pygmaeus",       "Pongo sp.",             "ambiguous"),
]
EXTRA = {"Cebuella pygmaea": "Callithrix_pygmaea", "Gorilla sp": "Gorilla_gorilla",
         "Tarsius sp": "Tarsius_syrichta", "Lagothrix lagothricha": "Lagothrix_lagotricha",
         "Callicebus moloch": "Plecturocebus_moloch"}

def sp_space(x): return re.sub(r'\s+', ' ', str(x).replace('_', ' ').strip())
def canon_species(x):
    s = sp_space(x)
    for tree, merge, kind in ALIAS:
        if kind != 'ambiguous' and s == merge: return tree
    if s in EXTRA: return EXTRA[s]
    return s.replace(' ', '_')
def to_merge_species(x):
    t = sp_space(x).replace(' ', '_')
    for tree, merge, kind in ALIAS:
        if t == tree: return merge
    return sp_space(t)

# ---------- read repo data ---------------------------------------------------
def read(path):
    with open(ROOT + path, newline='', encoding='utf-8-sig') as f:
        return list(csv.DictReader(f))

stephan = [r for r in read('data_raw/Stephan_primates.csv') if (r.get('Species') or '').strip()]
xwalk   = [r for r in read('metadata/qc_stephan/stephan_wide_crosswalk.csv')]
tips    = [x.strip() for x in re.findall(r'[(,]\s*([A-Za-z_][A-Za-z_. ]*)',
                                         open(ROOT + 'data_raw/species.nwk').read()) if x.strip()]
for r in stephan: r['tip'] = canon_species(r['Species'])

print('\n== 1. alias table ==')
check('canon_species is idempotent on all 59 tips',
      all(canon_species(t) == t for t in tips),
      [t for t in tips if canon_species(t) != t])
check('to_merge_species round-trips back to the tip (except ambiguous Pongo)',
      all(canon_species(to_merge_species(t)) == t for t in tips if not t.startswith('Pongo_')),
      [t for t in tips if not t.startswith('Pongo_') and canon_species(to_merge_species(t)) != t])
check('Cebuella_pygmaea folds onto the tip (the Matano-rename safety net)',
      canon_species('Cebuella_pygmaea') == 'Callithrix_pygmaea', canon_species('Cebuella_pygmaea'))
check("'Pongo sp.' is never auto-resolved to one tip",
      canon_species('Pongo sp.') == 'Pongo_sp.', canon_species('Pongo sp.'))
check('both Pongo tips map onto the single pooled merge row',
      to_merge_species('Pongo_abelii') == to_merge_species('Pongo_pygmaeus') == 'Pongo sp.')
check('space form and underscore form of a merge label agree',
      canon_species('Callicebus moloch') == canon_species('Callicebus_moloch') == 'Plecturocebus_moloch')

print('\n== 2. Stephan_primates.csv vs the tree ==')
check('59 rows', len(stephan) == 59, len(stephan))
check('Stephan tip set == species.nwk tip set exactly',
      set(r['tip'] for r in stephan) == set(tips),
      set(r['tip'] for r in stephan) ^ set(tips))
check('no duplicate tips after aliasing',
      len(set(r['tip'] for r in stephan)) == len(stephan))

print('\n== 3. crosswalk csv ==')
roles = Counter(r['role'] for r in xwalk)
check('65 pairs, 14 study3 + 1 predictor + 50 reference',
      len(xwalk) == 65 and roles['study3'] == 14 and roles['predictor'] == 1 and roles['reference'] == 50,
      dict(roles))
check('labels unique', len(set(r['label'] for r in xwalk)) == len(xwalk))
check('stephan_var unique', len(set(r['stephan_var'] for r in xwalk)) == len(xwalk))
sh = stephan[0].keys()
missing = [r['stephan_var'] for r in xwalk if r['role'] in ('study3', 'predictor') and r['stephan_var'] not in sh]
check('every compared stephan_var exists in Stephan_primates.csv', not missing, missing)
check('unit_scale parses', all(float(r['unit_scale']) > 0 for r in xwalk))
smap_pre = {r['tip']: r for r in stephan}

# ---------- build a merge-form mock of volumes_wide_select.csv ---------------
# The real wide file is on OneDrive (not mounted). Mock it by copying Stephan's
# own values across the crosswalk, in the MERGE namespace, so that:
#   * the alias layer is exercised for real (pooled "Pongo sp.", "Callicebus
#     moloch", "Lagothrix lagothricha", "Gorilla sp.", "Tarsius sp.")
#   * both sides hold IDENTICAL numbers by construction, so any status other
#     than "same" is a logic error and not a data disagreement -- except for the
#     pooled row, where one merge row must serve two tips with different values.
compared_x = [x for x in xwalk if x['role'] in ('study3', 'predictor')]
wide_cols = [x['wide_var'] for x in compared_x]
mock, pooled_loser = [], None
for t in tips:
    m = to_merge_species(t)
    if any(r['Species'] == m for r in mock):      # pooled: merge row already made
        pooled_loser = t                          # this tip must disagree with it
        continue
    row = {'Species': m}
    for x in compared_x:
        v = smap_pre[t].get(x['stephan_var'], '')
        row[x['wide_var']] = '' if str(v).strip() == '' else str(float(v) * float(x['unit_scale']))
    mock.append(row)
mock.append({'Species': 'Otolemur garnettii', **{c: '' for c in wide_cols}})  # extra, only-in-wide
for r in mock: r['merge_sp'] = sp_space(r['Species'])

print('\n== 4. species pairing (Stephan tip -> merge label -> wide row) ==')
wide_keys = set(r['merge_sp'] for r in mock)
pairs = [{'tip': r['tip'], 'merge_sp': to_merge_species(r['tip'])} for r in stephan]
pairs = [p for p in pairs if p['merge_sp'] in wide_keys]
cnt = Counter(p['merge_sp'] for p in pairs)
for p in pairs: p['pooled'] = cnt[p['merge_sp']] > 1
only_stephan = sorted(set(r['tip'] for r in stephan) - set(p['tip'] for p in pairs))
only_wide = sorted(k for k in wide_keys if k not in set(p['merge_sp'] for p in pairs))
check('all 59 Stephan rows pair', len(pairs) == 59, len(pairs))
check('nothing left only-in-Stephan', not only_stephan, only_stephan)
check('the 2 Pongo tips are flagged pooled, not silently collapsed',
      sorted(p['tip'] for p in pairs if p['pooled']) == ['Pongo_abelii', 'Pongo_pygmaeus'],
      [p['tip'] for p in pairs if p['pooled']])
check('a genuinely extra wide row is reported, not swallowed',
      only_wide == ['Otolemur garnettii'], only_wide)

print('\n== 5. regression: the damage the archived Matano script did ==')
broken = [dict(r, tip='Cebuella_pygmaea') if r['tip'] == 'Callithrix_pygmaea' else r for r in stephan]
# simulate the OLD behaviour: no alias layer, plain underscore matching
old = [r['tip'] for r in broken if sp_space(r['tip']) in wide_keys or r['tip'] in wide_keys]
check('WITHOUT the alias layer the renamed row is lost (this was the bug)',
      'Cebuella_pygmaea' not in old)
new = [canon_species(r['tip']) for r in broken]
new = [t for t in new if to_merge_species(t) in wide_keys]
check('WITH the alias layer the same file still pairs all 59', len(new) == 59, len(new))

print('\n== 6. cell comparison + anomaly detection ==')
TOL_SAME, TOL_FLAG = 1e-6, 0.01
def nm(x):
    try: return float(str(x).replace(',', ''))
    except Exception: return None
smap = {r['tip']: r for r in stephan}
wmap = {r['merge_sp']: r for r in mock}
cells = []
for x in xwalk:
    if x['role'] not in ('study3', 'predictor'): continue
    if x['stephan_var'] not in sh or x['wide_var'] not in wide_cols: continue
    for p in pairs:
        a, b = nm(smap[p['tip']].get(x['stephan_var'])), nm(wmap[p['merge_sp']].get(x['wide_var']))
        if a is not None: a *= float(x['unit_scale'])
        if a is None and b is None: continue
        rel = abs(a - b) / max(abs(a), abs(b), 1e-300) if (a is not None and b is not None) else None
        st = ('Stephan_only' if b is None else 'wide_only' if a is None
              else 'same' if rel <= TOL_SAME else 'rounding' if rel <= TOL_FLAG else 'DIFFERENT')
        ratio = b / a if (a and b) else None
        cells.append(dict(tip=p['tip'], label=x['label'], wide_var=x['wide_var'],
                          merge_sp=p['merge_sp'], a=a, b=b, rel=rel, status=st, ratio=ratio))
check('cells produced', len(cells) > 0, len(cells))
check('all 15 compared variables resolve on both sides',
      len(set(c['label'] for c in cells)) == 15,
      sorted(set(x['label'] for x in xwalk if x['role'] in ('study3', 'predictor'))
             - set(c['label'] for c in cells)))
st = Counter(c['status'] for c in cells)
print('        status counts:', dict(st))
# identical numbers by construction, so the ONLY legitimate disagreements are
# the second tip under the pooled "Pongo sp." row.
bad = sorted(set(c['tip'] for c in cells if c['status'] == 'DIFFERENT'))
check('no spurious DIFFERENT: every one is the pooled tip (%s)' % pooled_loser,
      bad in ([], [pooled_loser]), bad)
check('no spurious Stephan_only / wide_only outside the pooled tip',
      all(c['tip'] == pooled_loser for c in cells
          if c['status'] in ('Stephan_only', 'wide_only')),
      sorted(set(c['tip'] for c in cells if c['status'] in ('Stephan_only', 'wide_only')
                 and c['tip'] != pooled_loser)))
check('unit_scale is actually applied (cm3 -> mm3 pairs stay "same")',
      all(c['status'] in ('same', 'rounding') for c in cells
          if c['tip'] != pooled_loser and c['a'] and c['b']))

def near(x, t): return x is not None and abs(x / t - 1) < 0.05
def anomaly(r):
    if near(r, 2): return '~2x'
    if near(r, 0.5): return '~0.5x'
    if near(r, 1000): return '~1000x'
    if near(r, 0.001): return '~0.001x'
    return None
# inject known faults and confirm each is caught
inj = [(2.0, '~2x'), (0.5, '~0.5x'), (1000.0, '~1000x'), (0.001, '~0.001x'), (1.004, None)]
for factor, want in inj:
    got = anomaly(factor)
    check('ratio %-7s -> %-8s' % (factor, want), got == want, got)
c1 = next(c for c in cells if c['a'] and c['status'] == 'same')
c1b = dict(c1, b=c1['a'] * 2)
c1b['rel'] = abs(c1b['a'] - c1b['b']) / max(abs(c1b['a']), abs(c1b['b']))
check('a doubled cell is classified DIFFERENT and tagged ~2x',
      c1b['rel'] > TOL_FLAG and anomaly(c1b['b'] / c1b['a']) == '~2x')

print('\n== 7. wide-side source resolution ==')
def src_col(v, cols):
    for c in (v + '_Source', v + '_source', 'Source_' + v, re.sub(r'_Vol\.mm3$', '', v) + '_Source'):
        if c in cols: return c
    return None
check('no *_Source column in the snapshot -> falls back to the ledger',
      all(src_col(v, wide_cols) is None for v in set(c['wide_var'] for c in cells)))
withsrc = wide_cols + ['Amygdala_Vol.mm3_Source']
check('a *_Source column IS detected when present',
      src_col('Amygdala_Vol.mm3', withsrc) == 'Amygdala_Vol.mm3_Source')
check('the _Vol.mm3-stripped variant is detected too',
      src_col('Amygdala_Vol.mm3', wide_cols + ['Amygdala_Source']) == 'Amygdala_Source')

# synthetic ledger exercising the tie-break rule
ledger = [
    dict(Species='Homo sapiens', Variable='Amygdala_Vol.mm3', Value=1642.0, Source='Stephan_etal_1981_TableV',   Team='Selected_collection', Year=1981),
    dict(Species='Homo sapiens', Variable='Amygdala_Vol.mm3', Value=1642.0, Source='Bauernfeind_etal_2013',      Team='Selected_collection', Year=2013),
    dict(Species='Homo sapiens', Variable='Amygdala_Vol.mm3', Value=1642.0, Source='SomeOther_2020',             Team='Other_team',          Year=2020),
    dict(Species='Homo sapiens', Variable='Amygdala_Vol.mm3', Value=9999.0, Source='Wrong_value_1999',           Team='Selected_collection', Year=1999),
]
def resolve(species, var, value):
    j = [r for r in ledger if r['Species'] == species and r['Variable'] == var]
    hit = [r for r in j if abs(r['Value'] - value) <= 0.001 * max(abs(r['Value']), abs(value))]
    if not hit: return None
    sel = [r for r in hit if r['Team'] == 'Selected_collection']
    if sel: hit = sel
    hit.sort(key=lambda r: -r['Year'])
    return hit[0]['Source'] + (' (+%d tied)' % (len(hit) - 1) if len(hit) > 1 else '')
check('picks Selected_collection over a newer outside team, and reports the tie',
      resolve('Homo sapiens', 'Amygdala_Vol.mm3', 1642.0) == 'Bauernfeind_etal_2013 (+1 tied)',
      resolve('Homo sapiens', 'Amygdala_Vol.mm3', 1642.0))
check('a value with no ledger match yields no source rather than a wrong one',
      resolve('Homo sapiens', 'Amygdala_Vol.mm3', 123.0) is None)
check('ledger is keyed in the MERGE namespace, so tips must be converted first',
      to_merge_species('Plecturocebus_moloch') == 'Callicebus moloch')

print('\n== 8. Stephan-side source join ==')
refs = read('metadata/qc_stephan/Stephan_primates_references_long.csv')
rmap = {}
for r in refs:
    k = (canon_species(r['Species']), r['Stephan_column'])
    rmap.setdefault(k, r.get('preferred_reference') or r.get('status'))
check('references_long.csv species all canonicalise onto tips',
      set(canon_species(r['Species']) for r in refs) <= set(tips),
      sorted(set(canon_species(r['Species']) for r in refs) - set(tips)))
hit = sum(1 for c in cells if (c['tip'], next(x['stephan_var'] for x in xwalk if x['label'] == c['label'])) in rmap)
check('most compared cells get a Stephan-side source', hit / len(cells) > 0.5,
      '%d/%d (%.0f%%)' % (hit, len(cells), 100 * hit / len(cells)))
bycol = read('metadata/qc_stephan/Stephan_primates_references_by_column.csv')
asg = next(r for r in bycol if r['Stephan_column'] == 'ASG_Sousa')
check('ASG_Sousa really does resolve 0 cells to its stated source (the flagged column)',
      int(asg['n_resolved']) == 0, 'n_cells=%s resolved=%s differs=%s gap=%s'
      % (asg['n_cells'], asg['n_resolved'], asg['n_value_differs'], asg['n_gap']))

print('\n' + '=' * 60)
print('FAILED: %d' % len(FAIL))
for f in FAIL: print('   -', f)
sys.exit(1 if FAIL else 0)
