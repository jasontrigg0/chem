#!/usr/bin/env python
import re

import utils
import Queue
import collections
import copy

logger = utils.Logger()

class Molecule(object):
    def __init__(self):
        self.atoms = []
        self.atom2label = {} #for labels like [C:4]
        self.label2atom = {}

        #maintain two lists of bonds -- self.bonds for quick access
        #and self.ordered_bonds to store stereochemical information
        self.bonds = {} #atom1 -> atom2 -> bond_type (bonds are symmetric)
        self.ordered_bonds = {} #atom1 -> [(atom2,bond_type)]

        #when the ring_bonds are removed,
        #the rest of the nodes form a tree
        self.children = {}
        self.parent = {}
        self.ring_bonds = {} #{atom1 : {atom2 : label_number}}
    def copy(self):
        out = Molecule()
        out.atoms = copy.copy(self.atoms)
        out.atom2label = copy.copy(self.atom2label)
        out.label2atom = copy.copy(self.label2atom)

        #the self.bonds and self.ordered_bonds dictionaries
        #are nested, so copy.copy won't work
        #also deepcopy will copy the atom elements themselves, which
        #we don't want
        #so re-add the bonds:
        for atom,v in self.bonds.items():
            for atom2,bond in v.items():
                out.add_bond(atom,atom2,bond)

        out.children = copy.copy(self.children)
        out.parent = copy.copy(self.parent)
        out.ring_bonds = copy.copy(self.ring_bonds)
        return out
    def last_atom(self):
        if self.atoms:
            return self.atoms[-1]
        else:
            return None
    def add_bond(self, atom1, atom2, bond_type):
        if self.bonds.get(atom1,{}).get(atom2,None):
            return
        self.bonds.setdefault(atom1,{})[atom2] = bond_type
        self.bonds.setdefault(atom2,{})[atom1] = bond_type
        #adding bond
        self.ordered_bonds.setdefault(atom1,[]).append((atom2,bond_type))
        self.ordered_bonds.setdefault(atom2,[]).append((atom1,bond_type))
    def bond_summary(self, atom):
        for a2, bond_type in self.ordered_bonds.get(atom,[]):
            yield (str(a2), str(bond_type))
    def all_bonds(self, raw=False):
        out = []
        for a1 in self.atoms:
            for a2, bond_type in self.ordered_bonds.get(a1,[]):
                if raw:
                    out.append((a1, a2, bond_type))
                else:
                    out.append((str(a1), str(a2), str(bond_type)))
        return out
    def remove_atom(self, atom1):
        self.atoms.remove(atom1)
        self.bonds.pop(atom1)
        self.ordered_bonds.pop(atom1)
    def remove_bond(self, atom1, atom2, bond_type):
        if not self.bonds.get(atom1,{}).get(atom2,None) or\
           not self.bonds.get(atom2,{}).get(atom1,None):
            return
        self.bonds[atom1].pop(atom2)
        self.bonds[atom2].pop(atom1)
        self.ordered_bonds[atom1].remove((atom2, bond_type))
        self.ordered_bonds[atom2].remove((atom1, bond_type))
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
            for neighbor, bond_type in self.ordered_bonds[next_atom]:
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
            yield next_atom
    def __str__(self):
        return mol_to_smiles(self)

def map_labels(mol1, mol2):
    """use the labels in mol1, mol2 to
    generate a mapping from atoms in mol1 to atoms in mol2
    """
    dict_pairs = []
    for l in set(mol1.label2atom.keys() + mol2.label2atom.keys()):
        pair = (mol1.label2atom[l],mol2.label2atom[l])
        dict_pairs.append(pair)
    return dict(dict_pairs)
    

class Reaction(object):
    def __init__(self, smiles1, smiles2):
        self.rxn_reactants = smiles_to_mol(smiles1)
        self.rxn_products = smiles_to_mol(smiles2)
    def run(self, reactants_mol):
        mappings = find_substructure(self.rxn_reactants, reactants_mol)
        out = []
        for m in mappings:
            mol_to_transform = reactants_mol.copy()

            #remove all the reactant atom bonds:
            for a1 in self.rxn_reactants.atoms:
                for neighbor, bond_type in self.rxn_reactants.ordered_bonds[a1]:
                    mol_to_transform.remove_bond(m[a1], m[neighbor], bond_type)
                if not a1 in self.rxn_reactants.atom2label:
                    #TODO: delete it
                    mol_to_transform.remove_atom(m[a1])

            #generate mapping from rxn_products to reactants_mol using
            #the labels to map from rxn_products to rxn_reactants + 
            #the existing mapping from rxn_reactants to reactants_mol
            products_map = {}
            for a2 in self.rxn_products.atoms:
                if a2 in self.rxn_products.atom2label:
                    label = self.rxn_products.atom2label[a2]
                    reactant = self.rxn_reactants.label2atom[label]
                    products_map[a2] = m[reactant]
                else:
                    new_atom = Atom(str(a2))
                    mol_to_transform.atoms.append(new_atom)
                    products_map[a2] = new_atom

            #add all the reaction product bonds:
            for a2 in self.rxn_products.atoms:
                for neighbor, bond_type in self.rxn_products.ordered_bonds[a2]:
                    if not neighbor in products_map:
                        # print "missing from mapping!"
                        #create new atom and add to mapping
                        new_atom = Atom(str(neighbor))
                        mol_to_transform.atoms.append(new_atom)
                        products_map[neighbor] = new_atom
                    
                    mol_to_transform.add_bond(products_map[a2], products_map[neighbor], bond_type)
            out.append(mol_to_transform)
        return out
        
class Atom(object):
    """Name, Bonds/Neighbors, Orientation/Stereochemistry, Labels"""
    def __init__(self, name):
        self.name = name
    # def bond_summaries(self):
    #     return [(str(bond_type), str(a)) for bond_type, a in self.bonds]
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
            mol.label2atom[label] = next_atom
            mol.atom2label[next_atom] = label
        b = Bond(bond_string)
        if atom:
            mol.add_bond(atom, next_atom,b)
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
    parent = mol.parent.get(atom,None)
    out = ""
    if parent:
        out += str(mol.bonds[atom][parent])
    out += str(atom)
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

    def is_bond_subset(mol1, atom1, mol2, atom2):
        return is_subset(mol1.bond_summary(atom1), mol2.bond_summary(atom2))

    possibilities = {}
    for a1 in sub_mol.atoms:
        for a2 in mol.atoms:
            if str(a1) == str(a2) and is_bond_subset(sub_mol, a1, mol, a2):
                possibilities.setdefault(a1,[]).append(a2)

    for a1 in sub_mol.atoms:
        for a2 in mol.atoms:
            pass

    return match_substructure_helper(sub_mol, mol, possibilities, {}, {})
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
            for neighbor, bond_type in sub_mol.ordered_bonds[a1]:
                if not (matches[neighbor], bond_type) in mol.ordered_bonds[b1]:
                    return [] #no match!
        return [matches.copy()]

    #recursion
    out = []
    for a1 in sub_mol.atoms:
        if a1 in matches:
            continue
        if a1 not in possibilities:
            return []
        for m in possibilities[a1]:
            if m in matches.values(): #don't map two atoms from sub_mol to the same atom in mol
                continue
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

def test_substructure():
    m1 = smiles_to_mol("CCC(Br)CI")
    # m2 = smiles_to_mol('[C:1][C:2]([C:3])(Br)[C:4]') #>>[C:1][C:2]([C:3])=[C:4]')
    m2 = smiles_to_mol("C(Br)CI")
    assert(len(find_substructure(m2, m1)) == 1)

    m1 = smiles_to_mol("CCC(Br)C")
    # m2 = smiles_to_mol('[C:1][C:2]([C:3])(Br)[C:4]') #>>[C:1][C:2]([C:3])=[C:4]')
    m2 = smiles_to_mol("C(Br)C")
    assert(len(find_substructure(m2, m1)) == 2)
    # print find_substructure(m2, m1)

def test_reaction():
    m1 = smiles_to_mol("CC(C)(Br)CI")
    m2 = smiles_to_mol('[C:1][C:2]([C:3])(Br)[C:4]')
    m3 = smiles_to_mol('[C:1][C:2]([C:3])=[C:4]')
    rxn = Reaction('[C:1][C:2]([C:3])(Br)[C:4]','[C:1][C:2]([C:3])=[C:4]')
    print [str(out) for out in rxn.run(m1)]

def test_all():
    test_mol_to_smiles()
    test_substructure()
    
if __name__ == "__main__":
    test_reaction()
    # m1 = smiles_to_mol("CCC(Br)CI")
    # print mol_to_smiles(m1)
