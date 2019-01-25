#!/usr/bin/env python3
# coding=utf-8

import sys
import pybel
import math
import multiprocessing

sys.path.append('..')
sys.path.append('../../ms-tools')
from app.models import db
from mstools.formula import Formula

if sys.argv[1] == 'nist':
    from app.models import NistMolecule as Molecule, NistGroup as Group, NistMoleculeGroup as MoleculeGroup
elif sys.argv[1] == 'pubchem':
    from app.models_pubchem import PubchemMolecule as Molecule, PubchemGroup as Group, \
        PubchemMoleculeGroup as MoleculeGroup

smarts_bad = {
    'radicalC'                     : '[#6;v0,v1,v2,v3]',
    '*=*=*'                        : '*=*=*',
    '*#*~*#*'                      : '*#*~*#*',
    '[F,Cl,Br]~[!#6]'              : '[F,Cl,Br]~[!#6]',
    '*#*~[!#6]'                    : '*#*~[!#6]',
    '[NX2,NX4]'                    : '[NX2,NX4]',
    'O~N(~[!$([OX1])])~[!$([OX1])]': 'O~N(~[!$([OX1])])~[!$([OX1])]',
    'peroxide'                     : 'O~O',
    'N~N'                          : 'N~N',
    '[O,N]*[O,N;H1,H2]'            : '[O,N]*[O,N;H1,H2]',
    'C=C~[O,N;H1,H2]'              : 'C=C~[O,N;H1,H2]',
    'beta-dicarbonyl'              : 'O=C~*~C=O',
    'a=*'                          : 'a=*',
    'o'                            : 'o',
    '[n;r5]'                       : '[n;r5]',
    'pyridine-N-oxide'             : '[nX3;r6]',
    'triazine(zole)'               : '[$(nnn),$(nnan),$(nanan)]',
    '[R3]'                         : '[R3]',
    '[r3,r4;R2]'                   : '[r3,r4;R2]',
    '[r3,r4;#6X3]'                 : '[r3,r4;#6X3]',
    '[r3,r4]~[!#6]'                : '[r3,r4]~[!#6]',
    'nitrate'                      : 'O[NX3](~[OX1])~[OX1]',
    'amide'                        : 'O=C[NX3]',
    'acyl-halide'                  : 'O=C[F,Cl,Br]',
    'polybenzene'                  : 'c1ccc2c(c1)cccc2',
}
smarts_good = {
    'C=C~[O,N;H0]': 'C=C~[O,N;H0]',
    '[r5]~[!#6]'  : '[r5]~[!#6]',
    'ring3'       : '[r3]',
    'ring4'       : '[r4]',
    'ring5'       : '[r5]',
    'alkene'      : '[CX3]=[CX3]',
    'alkyne'      : '[CX2]#[CX2]',
    'benzene'     : 'c1ccccc1',
    # 'ether'       : 'COC',
    # 'Ph-ether'    : 'cO[C,c]',
    'ether'       : '[C,c]O[C,c]',
    'acetal'      : '[O;H0]C[O;H0]',
    'ketone'      : '[C,c]C(=O)[C,c]',
    'aldehyde'    : '[C,c]C(=O)[H]',
    # 'alcohol'     : 'C[OH]',
    # 'Ph-ol'       : 'c[OH]',
    'alcohol'     : '[C,c][OH]',
    'ester'       : 'C(=O)O[C,c;!$(C=O)]',
    'acid'        : 'C(=O)[OH]',
    # 'amine'       : 'C[NX3;!$(NC=O);!$(N~O);!$(Nc)]',
    # 'Ph-amine'    : 'c[NX3;!$(NC=O);!$(N~O)]',
    'amine'       : '[C,c][NX3;!$(NC=O);!$(N~O)]',
    'nitrile'     : '[NX1]#[CX2][C,c]',
    'nitro'       : '[C,c][NX3](~[OX1])~[OX1]',
    'N-aromatic'  : 'n',
    'halogen'     : '[F,Cl,Br]',
    # 'F'           : 'F',
    # 'Cl'          : 'Cl',
    # 'Br'          : 'Br',
    'CX4-X'       : '[CX4][F,Cl,Br]',
    'c-X'         : 'c[F,Cl,Br]',
    'C=C-X'       : 'C=C[F,Cl,Br]',
}
smarts_dict = dict(**smarts_bad, **smarts_good)

db.create_all()
db.session.query(MoleculeGroup).delete()
db.session.query(Group).delete()

for name, smarts in smarts_dict.items():
    group = Group(name=name, smarts=smarts)
    if name in smarts_bad.keys():
        group.bad = True
    db.session.add(group)

diol = Group(name='diol', smarts='diol')
hydrocarbon = Group(name='hydrocarbon', smarts='hydrocarbon')
alkane = Group(name='alkane', smarts='alkane')

db.session.add(diol)
db.session.add(hydrocarbon)
db.session.add(alkane)

db.session.commit()

smarts_id = dict()
for group in Group.query:
    smarts_id[group.smarts] = group.id


def match(smiles):
    group_ids = []
    if smiles.find('.') > -1:  # ignore complex
        return group_ids

    py_mol = pybel.readstring('smi', smiles)
    py_mol.removeh()
    formula = py_mol.formula
    atom_set = set(Formula.read(formula).atomdict.keys())
    if not atom_set <= {'C', 'H', 'O', 'N', 'F', 'Cl', 'Br'}:
        return group_ids

    if not ('C' in atom_set and {'H', 'F', 'Cl', 'Br'} & atom_set != set()):
        return group_ids

    if atom_set == {'C', 'H'}:
        group_ids.append(smarts_id['hydrocarbon'])

        for s in ['[CX2]', '[CX3]', 'c', '[#6;v0,v1,v2,v3]']:
            if pybel.Smarts(s).findall(py_mol) != []:
                break
        else:
            group_ids.append(smarts_id['alkane'])

    if pybel.Smarts('[OH]').findall(py_mol).__len__() > 1:
        group_ids.append(smarts_id['diol'])

    for name, smarts in smarts_dict.items():
        if pybel.Smarts(smarts).findall(py_mol):
            group_ids.append(smarts_id[smarts])

    return group_ids


def match_all(n_procs=8):
    print('Match group from smiles...')
    mols = db.session.query(Molecule).filter(Molecule.smiles != None).all()

    chunk_size = 10000
    for i in range(math.ceil(len(mols) / chunk_size)):
        sys.stdout.write('\r\t%i' % (i * chunk_size))

        mols_chunk = mols[i * chunk_size:(i + 1) * chunk_size]
        with multiprocessing.Pool(n_procs) as p:
            results = p.map(match, [mol.smiles for mol in mols_chunk])

        for mol, ids in zip(mols_chunk, results):
            mol_group_list = [MoleculeGroup(molecule_id=mol.id, group_id=id) for id in ids]
            db.session.add_all(mol_group_list)

    db.session.commit()


if __name__ == '__main__':
    match_all()
