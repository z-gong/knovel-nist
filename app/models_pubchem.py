from sqlalchemy.orm import relationship
from sqlalchemy import Column, Integer, Float, Text, String, ForeignKey, Boolean

from app import db


class PubchemGroup(db.Model):
    __bind_key__ = 'pubchem'
    __tablename__ = 'pubchem_group'
    id = Column(Integer, primary_key=True)
    name = Column(String(100))
    smarts = Column(String(100))
    bad = Column(Boolean, default=False)

    molecules = relationship('PubchemMolecule', secondary='pubchem_molecule_group')

    def __repr__(self):
        return '<Group: %s %s>' % (self.name, self.smarts)


class PubchemMolecule(db.Model):
    __bind_key__ = 'pubchem'
    __tablename__ = 'pubchem_molecule'
    id = Column(Integer, primary_key=True)
    cid = Column(Integer, unique=True)
    cas = Column(String(100))  # CAS could be duplicated for isomers
    formula = Column(String(100))
    name = Column(Text)
    smiles = Column(Text)  # Convert to canonical SMILES
    inchi = Column(Text)  # Original InChI from PubChem SDF. Leave it empty to reduce the size of db
    weight = Column(Float)
    n_heavy = Column(Integer)
    remark = Column(Text)

    groups = relationship('PubchemGroup', secondary='pubchem_molecule_group')

    def __repr__(self):
        return '<PubchemMolecule: %s %s %s>' % (self.formula, self.smiles, self.name)

    def __str__(self):
        return '%i\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
            self.id, self.formula, 'None', self.smiles, 'None', 'None', 'None', self.cid)


class PubchemMoleculeGroup(db.Model):
    __bind_key__ = 'pubchem'
    __tablename__ = 'pubchem_molecule_group'
    molecule_id = Column(Integer, ForeignKey(PubchemMolecule.id), primary_key=True)
    group_id = Column(Integer, ForeignKey(PubchemGroup.id), primary_key=True)
