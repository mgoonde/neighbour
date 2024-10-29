# Description

Example program which reads an xyz structure, computes the neighbour list with `rcut=6.0`,
and outputs the `list` and `veclist` for some atom index `idx=3`.

# Compile

First compile the `m_neighbour.f90` in the root directory (```cd ../ && sh compile.sh && cd -```), then:

```bash
sh compile.sh
```

# Run

```bash
./main.x < pt_eam.xyz
```
