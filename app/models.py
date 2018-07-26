""" Database engine
"""

from scipy import interpolate
import json
from sqlalchemy.orm import relationship
from sqlalchemy import Column, Integer, Float, Text, String, ForeignKey, Boolean

from . import db


class NistGroup(db.Model):
    __bind_key__ = 'nist'
    __tablename__ = 'nist_group'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    smarts = Column(String(100))
    bad = Column(Boolean, default=False)

    molecules = relationship('NistMolecule', secondary='nist_molecule_group')

    def __repr__(self):
        return '<Group: %s %s>' % (self.name, self.smarts)


class NistMolecule(db.Model):
    __bind_key__ = 'nist'
    __tablename__ = 'nist_molecule'
    id = Column(Integer, primary_key=True)
    content_id = Column(String(100), unique=True)
    cas = Column(String(100))  # CAS could be duplicated for isomers
    formula = Column(String(100))
    name = Column(Text)
    smiles = Column(Text)  # convert to canonical SMILES
    inchi = Column(Text)  # original InChI from knovel-nist website
    weight = Column(Float)
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
    remark = Column(Text)
    n_heavy = Column(Integer)
    constant_inserted = Column(Boolean, default=False)
    data_inserted = Column(Boolean, default=False)
    spline_inserted = Column(Boolean, default=False)

    datas = relationship('NistData', lazy='dynamic')
    splines = relationship('NistSpline', lazy='dynamic')

    groups = relationship('NistGroup', secondary='nist_molecule_group')

    def __repr__(self):
        return '<NistMolecule: %s %s %s>' % (self.formula, self.smiles, self.name)

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
    def Tm25(self):
        if self.tt == None:
            return 300
        else:
            return max(self.tt + 25, 200)

    @property
    def Tcx8(self):
        if self.tc == None:
            return None
        else:
            return self.tc * 0.8

    @property
    def groups_str(self):
        groups = [group.name for group in self.groups]
        return ', '.join(groups)

    def get_Pvap(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'pvap-lg').first()
        if spline == None:
            return None
        pvap = spline.get_data(T)[0]
        if pvap != None:
            pvap /= 100
        return pvap  # bar

    def get_density(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'density-lg').first()
        if spline == None:
            return None
        dens = spline.get_data(T)[0]
        if dens != None:
            dens /= 1000
        return dens  # g/mL

    def get_Hvap(self, T):
        spline = self.splines.join(NistProperty).filter(NistProperty.name == 'hvap-lg').first()
        if spline == None:
            return None
        hvap = spline.get_data(T)[0]
        return hvap  # kJ/mol

    def update_n_heavy(self):
        import pybel
        m = pybel.readstring('smi', self.smiles)
        self.n_heavy = m.OBMol.NumHvyAtoms()


class NistMoleculeGroup(db.Model):
    __bind_key__ = 'nist'
    __tablename__ = 'nist_molecule_group'
    molecule_id = Column(Integer, ForeignKey(NistMolecule.id), primary_key=True)
    group_id = Column(Integer, ForeignKey(NistGroup.id), primary_key=True)


class NistProperty(db.Model):
    __bind_key__ = 'nist'
    __tablename__ = 'nist_property'
    id = Column(Integer, primary_key=True)
    name = Column(Text)
    property_id = Column(String(100))
    phase_id = Column(String(100))

    datas = relationship('NistData', lazy='dynamic')
    splines = relationship('NistSpline', lazy='dynamic')


class NistData(db.Model):
    __bind_key__ = 'nist'
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
    __bind_key__ = 'nist'
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
