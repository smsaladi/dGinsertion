## Regression test datasets

By using proteins in the *E. coli* K-12 proteome that were identified to have
their C-terminus in the cytoplasm yields a good set of TMs and loops to
calculate dG_insertion values for. `residuelocation.py` is from
ml-expression/v1.0

```bash
python ~/ml-expression/scripts/residuelocation.py Daley_gfp.fna.faa Daley_gfp.fna.phobiuslong Daley_gfp.fna.mps Daley_gfp.fna.sps

cat Daley_gfp.fna.mps Daley_gfp.fna.sps | grep -v '>' | tr '_' '\n' | tr -d '*' | grep -v '^$' > Daley_gfp.segs

cat Daley_gfp.segs | ../calc_dG_v2.pl > Daley_gfp.segs.dGpred_local
```

Because changes were made to the local code without formal testing, the output
from both the local code as well as the online [service](http://dgpred.cbr.su.se)
are saved.
