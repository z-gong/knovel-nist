import math
from collections import OrderedDict
from flask import render_template, request, redirect, url_for
from . import main
from ..models import *
from sqlalchemy import or_
import pybel


@main.route('/')
def index():
    groups = NistGroup.query
    return render_template('index.html', groups=groups)


@main.route('/group')
def show_molecules_in_groups():
    group_least = int(request.args.get('group_least'))
    group_most = int(request.args.get('group_most'))
    heavy_least = int(request.args.get('heavy_least'))
    heavy_most = int(request.args.get('heavy_most'))
    include = request.args.get('include')
    exclude = request.args.get('exclude')

    molecules = NistMolecule.query.filter(NistMolecule.smiles != None) \
        .filter(NistMolecule.n_heavy >= heavy_least) \
        .filter(NistMolecule.n_heavy <= heavy_most)

    if include is not None and include != '':
        group_ids = [int(i) for i in include.split(',')]
        for gid in group_ids:
            molecules = molecules.filter(NistMolecule.groups.any(NistGroup.id == gid))

    mol_list = molecules.all()
    if exclude is not None and exclude != '':
        exclude_ids = [int(i) for i in exclude.split(',')]
        for gid in exclude_ids:
            exclude_molecules = molecules.filter(NistMolecule.groups.any(NistGroup.id == gid)).all()
            mol_list = [mol for mol in mol_list if mol not in exclude_molecules]

    mol_list = filter(lambda x: len(x.groups) >= group_least, mol_list)
    mol_list = filter(lambda x: len(x.groups) <= group_most, mol_list)

    mol_list = sorted(mol_list, key=lambda x: x.n_heavy)

    return render_template('molecules.html', molecules=mol_list, T=298)


@main.route('/search/<mol_list_str>')
def search(mol_list_str):
    mol_list = mol_list_str.strip().split(',')
    for i, molstr in enumerate(mol_list):
        try:
            mol_list[i] = pybel.readstring('smi', molstr).write('can').strip()
        except:
            pass

    molecules = NistMolecule.query.filter(
        or_(
            NistMolecule.smiles.in_(mol_list),
            NistMolecule.formula.in_(mol_list),
            NistMolecule.name.in_(mol_list),
            NistMolecule.cas.in_(mol_list),
        )).all()
    return render_template('molecules.html', molecules=molecules, T=298)


@main.route('/train/<molecules>')
def train(molecules):
    mol_ids = list(map(int, molecules.split(',')))
    molecules = NistMolecule.query.filter(NistMolecule.id.in_(mol_ids)).all()
    return render_template('train.html', molecules=sorted(molecules, key=lambda x: x.n_heavy))


@main.route('/msdserver/<molecules>')
def msdserver(molecules):
    mol_ids = list(map(int, molecules.split(',')))
    molecules = []

    # sqlite has limit on the IN expression
    for i in range(math.ceil(len(mol_ids) / 500)):
        mols = NistMolecule.query.filter(NistMolecule.id.in_(mol_ids[i * 500:(i + 1) * 500]))
        Tvap = request.args.get('Tvap')
        if Tvap is not None:
            mols = mols.filter(NistMolecule.tb < int(Tvap) + 0.5)
        molecules += mols.all()

    molecules = sorted(molecules, key=lambda x: x.n_heavy)

    return render_template('msdserver.html', molecules=molecules)


@main.route('/validate/<molecules>')
def validate(molecules):
    mol_ids = list(map(int, molecules.split(',')))
    molecules = NistMolecule.query.filter(NistMolecule.id.in_(mol_ids))

    mol_T_dict = OrderedDict()
    for mol in molecules:
        if mol.Tfus is not None:
            Tmin = mol.Tfus
        else:
            Tmin = 200
        if mol.Tc is not None:
            Tmax = mol.Tc
        elif mol.Tvap is not None:
            Tmax = mol.Tvap + 100
        else:
            Tmax = 600

        dT = 25
        Tmin = math.ceil(Tmin / dT) * dT
        Tmax = math.floor(Tmax / dT) * dT
        T_list = range(Tmin, Tmax + 1, dT)
        mol_T_dict[mol] = T_list

    return render_template('validate.html', mol_T_dict=mol_T_dict)
