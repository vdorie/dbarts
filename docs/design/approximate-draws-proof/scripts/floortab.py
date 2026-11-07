"""LaTeX rows of Table tab:floor from floor_z0.out (one row per (b, psi))."""
import sys
def e(v):
    m, x = f'{float(v):.1e}'.split('e')
    return f'${m}\\mathrm{{e}}{{{int(x)}}}$'
print(r'\begin{tabular}{rrrrrrrrrr}\toprule')
print(r'& & & \multicolumn{2}{c}{untilted, 53 bits} & \multicolumn{4}{c}{tilted} & \\'
      r'\cmidrule(lr){4-5}\cmidrule(lr){6-9}')
print(r'$b$ & $|\psi|$ & $K$ & $\log\fz$ & LR & $\log f$ & $\log q$ & LR & $\log f$, 106 bits & $x$\\\midrule')
for line in open(sys.argv[1]):
    kv = dict(t.split('=', 1) for t in line.split() if '=' in t)
    if 'f53' not in kv: continue
    q = kv.get('q53L'); lr = kv.get('LR53L'); lro = kv.get('LR_old')
    print(f"{float(kv['b']):g} & {float(kv['psi']):g} & {kv.get('K','')} & {e(kv['f_old'])} & "
          f"{e(lro) if lro else '--'} & {e(kv['f53L'])} & {e(q) if q else '--'} & {e(lr) if lr else '--'} & "
          f"{e(kv['f106'])} & {float(kv['x']):.4g}\\\\")
print(r'\bottomrule\end{tabular}')
