import json
from scipy import interpolate
from sqlalchemy.orm import relationship
from sqlalchemy import Column, Integer, Float, Text, String, ForeignKey, Boolean
from openbabel import pybel
from .formula import Formula
from . import db


def ri(f):
    '''
    round float to int.
    None will still be None.
    '''
    if f is None:
        return None
    return int(round(f))


class NistGroup(db.Model):
    __tablename__ = 'nist_group'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    smarts = Column(String(100))
    bad = Column(Boolean, default=False)

    molecules = relationship('NistMolecule', secondary='nist_molecule_group')

    def __repr__(self):
        return '<Group: %s %s>' % (self.name, self.smarts)


class NistMolecule(db.Model):
    __tablename__ = 'nist_molecule'
    id = Column(Integer, primary_key=True)
    content_id = Column(String(100), unique=True)
    cas = Column(String(100))  # CAS could be duplicated for isomers or could be None
    formula = Column(String(100))
    name = Column(Text)
    smiles = Column(Text)  # convert to canonical SMILES
    inchi = Column(Text)  # original InChI from knovel-nist website
    weight = Column(Float)
    n_heavy = Column(Integer)
    n_nothx = Column(Integer)
    remark = Column(Text)

    tc = Column(Float)  # K
    pc = Column(Float)  # kPa
    dc = Column(Float)  # kg/m^3
    tb = Column(Float)  # K Boiling point
    tt = Column(Float)  # K Triple point
    hfus = Column(Float)  # kJ/mol
    tc_u = Column(Float)  # K
    pc_u = Column(Float)  # kPa
    dc_u = Column(Float)  # kg/m^3
    tb_u = Column(Float)  # K
    tt_u = Column(Float)  # K
    hfus_u = Column(Float)  # kJ/mol
    tc_has_exp = Column(Boolean, default=False)
    pc_has_exp = Column(Boolean, default=False)
    dc_has_exp = Column(Boolean, default=False)
    tb_has_exp = Column(Boolean, default=False)
    tt_has_exp = Column(Boolean, default=False)
    hfus_has_exp = Column(Boolean, default=False)

    constant_inserted = Column(Boolean, default=False)
    data_inserted = Column(Boolean, default=False)
    spline_inserted = Column(Boolean, default=False)

    datas = relationship('NistData', lazy='dynamic')
    splines = relationship('NistSpline', lazy='dynamic')

    groups = relationship('NistGroup', secondary='nist_molecule_group')

    def __repr__(self):
        return '<NistMolecule: %s %s %s>' % (self.formula, self.smiles, self.name)

    def __str__(self):
        return '%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (self.id, self.formula, self.cas or 'None', self.smiles,
                                                   ri(self.tt) or 'None', ri(self.tb) or 'None', ri(self.tc) or 'None',
                                                   self.content_id)

    @property
    def Tfus(self):
        return self.tt

    @property
    def Tvap(self):
        return self.tb

    @property
    def Tc(self):
        return self.tc

    @property
    def Pc(self):
        if self.pc is None:
            return None

        return self.pc / 100  # bar

    @property
    def Dc(self):
        if self.dc is None:
            return None

        return self.dc / 1000  # g/mL

    @property
    def Tm25(self):
        if self.tt is None:
            if self.tb is None:
                return 298
            return self.tb * 0.4 + 100
        else:
            return max(self.tt + 25, 200)

    @property
    def Tcx8(self):
        if self.tc is None:
            return None
        else:
            return self.tc * 0.8

    @property
    def groups_str(self):
        groups = [group.name for group in self.groups]
        return ', '.join(groups)

    def get_pvap(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'pvap-lg').first()
        if spline is None:
            return None
        pvap = spline.get_data(T)[0]
        if pvap is not None:
            pvap /= 100
        return pvap  # bar

    def get_density(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'density-lg').first()
        if spline is None:
            return None
        dens = spline.get_data(T)[0]
        if dens is not None:
            dens /= 1000
        return dens  # g/mL

    def get_hvap(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'hvap-lg').first()
        if spline is None:
            return None
        hvap = spline.get_data(T)[0]
        return hvap  # kJ/mol

    def get_cp(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'cp-lg').first()
        if spline is None:
            return None
        cp = spline.get_data(T)[0]
        return cp  # J/mol.K

    def get_viscosity(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'viscosity-lg').first()
        if spline is None:
            return None
        vis = spline.get_data(T)[0]
        if vis is not None:
            vis *= 1000
        return vis  # mPa.s

    def get_st(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'st-lg').first()
        if spline is None:
            return None
        st = spline.get_data(T)[0]
        if st is not None:
            st *= 1000
        return st  # mN/m

    def update_n_heavy(self):
        m = pybel.readstring('smi', self.smiles)
        self.n_heavy = m.OBMol.NumHvyAtoms()

        f = Formula(self.formula)
        N = 0
        for a, n in f.atoms.items():
            if a not in ('H', 'F', 'Cl', 'Br'):
                N += n
        self.n_nothx = N


class NistMoleculeGroup(db.Model):
    __tablename__ = 'nist_molecule_group'
    molecule_id = Column(Integer, ForeignKey(NistMolecule.id), primary_key=True)
    group_id = Column(Integer, ForeignKey(NistGroup.id), primary_key=True)


class NistProperty(db.Model):
    __tablename__ = 'nist_property'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    property_id = Column(String(100))
    phase_id = Column(String(100))

    datas = relationship('NistData', lazy='dynamic')
    splines = relationship('NistSpline', lazy='dynamic')


class NistHasData(db.Model):
    __tablename__ = 'nist_has_data'
    id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(NistMolecule.id))
    property_id = Column(Integer, ForeignKey(NistProperty.id))
    property_name = Column(String(100))  # save property name for convenience
    has_exp = Column(Boolean, default=False)  # has experimental data
    has_rec = Column(Boolean, default=False)  # has recommended data

    molecule = relationship(NistMolecule)
    property = relationship(NistProperty)


class NistData(db.Model):
    __tablename__ = 'nist_data'
    id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(NistMolecule.id))
    property_id = Column(Integer, ForeignKey(NistProperty.id))
    t = Column(Float)  # K
    p = Column(Float)  # kPa
    value = Column(Float)
    uncertainty = Column(Float)

    molecule = relationship(NistMolecule)
    property = relationship(NistProperty)


class NistSpline(db.Model):
    __tablename__ = 'nist_spline'
    id = Column(Integer, primary_key=True)
    molecule_id = Column(Integer, ForeignKey(NistMolecule.id))
    property_id = Column(Integer, ForeignKey(NistProperty.id))
    t_min = Column(Float)
    t_max = Column(Float)
    coef_v = Column(Text)  # value
    coef_u = Column(Text)  # uncertainty

    molecule = relationship(NistMolecule)
    property = relationship(NistProperty)

    def get_data(self, T):
        if T < self.t_min or T > self.t_max:
            return None, None

        coef_v = json.loads(self.coef_v)
        coef_u = json.loads(self.coef_u)

        v = interpolate.splev(T, coef_v)
        u = interpolate.splev(T, coef_u)
        return float(v), float(u)
