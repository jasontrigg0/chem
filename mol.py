#!/usr/bin/env python
import re

import utils
import Queue
import collections

logger = utils.Logger()

class Molecule(object):
    def __init__(self):
        self.atoms = []
        
        #when the ring_bonds are removed,
        #the rest of the nodes form a tree
        self.children = {}
        self.parent = {}
        self.ring_bonds = {} #{atom1 : {atom2 : label_number}}
    def last_atom(self):
        if self.atoms:
            return self.atoms[-1]
        else:
            return None
    def setup_dfs(self):
        self.children = {}
        self.parent = {}
        self.ring_bonds = {}
        if not self.atoms:
            return
        done = set()
        stack = collections.deque()
        stack.append(self.atoms[0])
        while stack:
            next_atom = stack.pop()
            for bond_type, neighbor in next_atom.bonds[::-1]:
                if neighbor in done and next_atom not in self.children.get(neighbor,[]):
                    self.ring_bonds.setdefault(next_atom,{}).setdefault(neighbor,{})
                    continue
                if neighbor in done and next_atom in self.children.get(neighbor,[]):
                    self.parent[next_atom] = neighbor
                else:
                    self.children.setdefault(next_atom,[]).append(neighbor)
                if neighbor not in done:
                    stack.append(neighbor)
            done.add(next_atom)
    def dfs(self, start_node=None):
        if not start_node and self.atoms:
            start_node = self.atoms[0]
        else:
            return
        stack = collections.deque()
        stack.append(start_node)
        while stack:
            next_atom = stack.pop()
            for a in self.children.get(next_atom,[])[::-1]:
                stack.append(a)
                print "appending... " + str(a)
            yield next_atom
    def __str__(self):
        """Convert to SMILES"""
        return str([str(a) for a in self.atoms])
    
class Atom(object):
    """Name, Bonds/Neighbors, Orientation/Stereochemistry, Labels"""
    def __init__(self, name):
        self.name = name
        self.bonds = [] #[(bond_type, Atom)]
    def add_bond(self, atom, bond):
        self.bonds.append((atom,bond))
    def __str__(self):
        return self.name

def test():
    pass

def parse_smiles(smiles):
    pass

class Bond(object):
    def __init__(self, bond_string):
        self.bond_string = bond_string
        if not bond_string:
            self.bond_cnt = 1
        elif bond_string == "=":
            self.bond_cnt = 2
        elif bond_string == "#":
            self.bond_cnt = 3
        else:
            raise
    def __str__(self):
        return str(self.bond_string)

    
def take_next(smiles):
    if smiles[0] == "(":
        #TODO: below regex will fail on nested branches
        branch, rest = re.findall(r"^(\(.*?\))(.*)$",smiles)[0]
        return branch, rest

    re_MAIN = r"^([{BOND_CHARS}]?(Cl|Br|O|I|C)\d*)(.*)$".format(**globals())
    m_main = re.findall(re_MAIN,smiles)
    if m_main:
        main, _, rest = m_main[0]
        return main, rest

    re_BRACKET = r"^([{BOND_CHARS}]?\[.*?\])(.*)$".format(**globals())
    m_bracket = re.findall(re_BRACKET, smiles)
    if m_bracket:
        bracket, rest = m_bracket[0]
        return bracket, rest

    # for char in smiles_list:
    #     if char not in bond_chars:
    #         next_atom = Atom(char)
    #         if last_atom:
    #             last_atom.add_bond(next_atom, BONDS.SINGLE)
    #             next_atom.add_bond(last_atom, BONDS.SINGLE)
    #         last_atom = next_atom

def get_chunks(smiles):
    while smiles:
        chunk, smiles = take_next(smiles)
        if not chunk: raise
        yield chunk

BOND_CHARS = "~=#"

def smiles_to_mol(smiles):
    mol = Molecule()
    mol = append_smiles_clause(smiles, mol)
    return mol

def append_smiles_clause(smiles_clause, mol, atom = None):
    if atom == None:
        atom = mol.last_atom()
    for chunk in get_chunks(smiles_clause):
        print "Current smiles: " + mol_to_smiles(mol)
        print "chunk: " + chunk
        if chunk[0] == "(" and chunk[-1] == ")":
            #append branch
            clause = chunk[1:-1]
            mol = append_smiles_clause(clause, mol, atom)
            continue
        
        if "[" in chunk:
            bond_string, main = re.findall(r"([{BOND_CHARS}]?)\[(.*)\]".format(**globals()),chunk)[0]
            #add
        else:
            bond_string, main = re.findall(r"([{BOND_CHARS}]?)(.*)".format(**globals()),chunk)[0]
        print "bond: " + bond_string
        print "main: " + main
        element, rest = parse_element(main)
        next_atom = Atom(element)
        b = Bond(bond_string)
        if atom:
            atom.bonds.append([b, next_atom])
            next_atom.bonds.append([b, atom])
        mol.atoms.append(next_atom)
        atom = mol.last_atom()
    return mol


def parse_element(main):
    atom_types = ["Cl","C","O","Br","I"]
    for a in atom_types:
        if main.startswith(a):
            return a, main[len(a):]
    return None
        
def add_atom(mol, next_atom):
    last_atom = mol.last_atom

def mol_to_smiles(mol):
    mol.setup_dfs()
    if not mol.atoms:
        return ""
    return smiles_helper(mol,mol.atoms[0])

def smiles_helper(mol, atom):
    """Compute smiles of a subsection of the molecule rooted at 'atom'"""
    children = mol.children.get(atom,[])
    out = str(atom)
    if len(children) > 1:
        for c in children[:-1]:
            out += "(" + smiles_helper(mol,c) + ")"
    if len(children) > 0:
        out += smiles_helper(mol,children[-1])
    return out

def test_chunks():
    s = "CCCC=O"
    s = '[CH3][CH2]C#[CH]'
    print s
    while s:
        chunk, s = take_next(s)
        print chunk
    
if __name__ == "__main__":
    s = '[CH3][CH2]C#[CH]'
    s = 'CC(Br)CC'
    m = smiles_to_mol(s)
    m.setup_dfs()
    print m
    print mol_to_smiles(m)
