#!/usr/bin/env python
import re

import utils
import Queue
import collections

logger = utils.Logger()

class Molecule(object):
    def __init__(self):
        self.atoms = []
        self.labels = {} #number -> atom
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
            for bond_type, neighbor in next_atom.bonds:
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
            for a in self.children.get(next_atom,[]):
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
    def bond_summaries(self):
        return [(str(bond_type), str(a)) for bond_type, a in self.bonds]
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
    def __eq__(self, other):
        return str(self.bond_string) == str(other.bond_string)

    
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
        element, rest = parse_element(main)
        next_atom = Atom(element)
        if ":" in rest: #C:1
            label = rest.split(":")[1]
            mol.labels[label] = next_atom
        b = Bond(bond_string)
        if atom:
            atom.bonds.append((b, next_atom))
            next_atom.bonds.append((b, atom))
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

def smarts2reaction(smarts, mol):
    reactants, products = smarts.split(">>",1)
    for chunk in get_chunks(reactants):
        print chunk

def counter(l, key=lambda x: x):
    cnts = {}
    for i in l:
        k = key(i)
        cnts[k] = cnts.get(k,0) + 1
    return cnts

def is_subset(list1, list2):
    return is_counter_subset(counter(list1),counter(list2))

def is_counter_subset(counter1, counter2):
    #check that counter1 < counter2
    for k1,v1 in counter1.items():
        if v1 > counter2.get(k1,0):
            return False
    return True

def find_substructure(sub_mol, mol):
    #atom -> (bond_type, atom_type) -> cnt
    # mol_cnts = dict((a,counter(a.bond_summaries())) for a in mol.atoms)

    def is_bond_subset(atom1, atom2):
        return is_subset(atom1.bond_summaries(), atom2.bond_summaries())

    possibilities = {}
    for a1 in sub_mol.atoms:
        for a2 in mol.atoms:
            if str(a1) == str(a2) and is_bond_subset(a1, a2):
                print a1.bond_summaries(), a2.bond_summaries()
                possibilities.setdefault(a1,[]).append(a2)

    for a1 in sub_mol.atoms:
        for a2 in mol.atoms:
            pass

    print match_substructure_helper(sub_mol, mol, possibilities, {}, {})
    #match the first atom
    # q = collections.deque()
    # done = {}
    # q.append(matches.keys()[0])
    # while q:
    #     next_atom = q.popleft()
    #     if next_atom in done:
    #         continue
    #     #try to match this atom
    #     for bond_type, neighbor in next_atom.bonds:
    #         pass
    #first check which atom from sub_mol is the hardest to match in mol 
    # for a in sub_mol.atoms:
    #     cnts = counter(a.bond_summaries())
    #     matches = []
    #     for b in mol_cnts:
    #         for k,v in cnts.items():
    # pass

def match_substructure_helper(sub_mol, mol, possibilities, matches, done):
    """matches = [(atom from sub_mol, mapped atom from mol)]"""
    #base case
    #all done matching -- check that it works!
    if len(matches) == len(sub_mol.atoms):
        for a1 in sub_mol.atoms:
            b1 = matches[a1]
            for bond_type, neighbor in a1.bonds:
                if not (bond_type, matches[neighbor]) in b1.bonds:
                    return [] #no match!
        return [matches.copy()]

    #recursion
    out = []
    for a1 in sub_mol.atoms:
        if a1 in matches:
            continue
        for m in possibilities[a1]:
            matches[a1] = m
            out += match_substructure_helper(sub_mol, mol, possibilities, matches, done)
            matches.pop(a1)
        break #only process one atom from sub_mol.atoms, then break
    return out

def run_reaction(react_smiles, product_smiles, mol):
    pass
    
def test_chunks():
    s = "CCCC=O"
    s = '[CH3][CH2]C#[CH]'
    print s
    while s:
        chunk, s = take_next(s)
        print chunk

def test_mol_to_smiles():
    s = '[CH3][CH2]C#[CH]'
    s = 'CC(Br)CC'
    m = smiles_to_mol(s)
    s_out = mol_to_smiles(m)
    assert(s_out == s)

if __name__ == "__main__":
    m1 = smiles_to_mol("CCC(Br)CI")
    # m2 = smiles_to_mol('[C:1][C:2]([C:3])(Br)[C:4]') #>>[C:1][C:2]([C:3])=[C:4]')
    m2 = smiles_to_mol("C(Br)CI")
    find_substructure(m2, m1)
