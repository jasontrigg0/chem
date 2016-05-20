#!/usr/bin/env python
import re

import utils
import Queue
import collections
import copy
import itertools

from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem

logger = utils.Logger()

class Structure(object):
    """Representation of a molecular structure:
    """
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
        
    def append_smiles_clause(self, smiles_clause, atom = None, template=False):
        if atom == None:
            atom = self.last_atom()
        for chunk in get_chunks(smiles_clause):
            if chunk[0] == "(" and chunk[-1] == ")":
                #append branch
                clause = chunk[1:-1]
                self.append_smiles_clause(clause, atom, template=template)
                continue
            if "[" in chunk:
                bond_string, main = re.findall(r"([{BOND_CHARS}]?)\[(.*)\]".format(**globals()),chunk)[0]
                #add
            else:
                bond_string, main = re.findall(r"([{BOND_CHARS}]?)(.*)".format(**globals()),chunk)[0]
            atom_string, label = parse_atom(main)
            if template:
                next_atom = AtomPattern(atom_string, self)
            else:
                next_atom = Atom(atom_string, self)
            if label: #C:1
                self.label2atom[label] = next_atom
                self.atom2label[next_atom] = label
            b = Bond(bond_string)
            if atom:
                self.add_bond(atom, next_atom,b)
            self.atoms.append(next_atom)
            atom = self.last_atom()
        
    def copy(self):
        out = Structure()
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
        out = []
        for a2, bond_type in self.ordered_bonds.get(atom,[]):
            out.append((a2, bond_type))
        return out
    def hydrogen_cnt(self, atom):
        pass
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
    
class Molecule(Structure):
    """
    Molecule has all the features of a Structure along with implied Hydrogen atoms
    """
    def __init__(self,smiles):
        super(Molecule,self).__init__()
        #parse the smiles clause and add atom-by-atom
        self.append_smiles_clause(smiles, template=False)
    
class Substructure(Structure):
    """
    """
    def __init__(self,smiles):
        super(Substructure,self).__init__()
        self.append_smiles_clause(smiles, template=True)

        
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
    def __init__(self, smarts):
        smiles1, smiles2 = smarts.split(">>",1)
        self.reactant_template = Substructure(smiles1)
        self.product_template = Substructure(smiles2)
    def run(self, reactants_tuple):
        reactants_mol = reactants_tuple[0]
        mappings = find_substructure(self.reactant_template, reactants_mol)
        out = []
        for m in mappings:
            mol_to_transform = reactants_mol.copy()

            #remove all the reactant atom bonds:
            for a1 in self.reactant_template.atoms:
                for neighbor, bond_type in self.reactant_template.ordered_bonds[a1]:
                    mol_to_transform.remove_bond(m[a1], m[neighbor], bond_type)
                if not a1 in self.reactant_template.atom2label:
                    #TODO: delete it
                    mol_to_transform.remove_atom(m[a1])

            #generate mapping from product_template to reactants_mol using
            #the labels to map from product_template to reactant_template + 
            #the existing mapping from reactant_template to reactants_mol
            products_map = {}
            for a2 in self.product_template.atoms:
                if a2 in self.product_template.atom2label:
                    label = self.product_template.atom2label[a2]
                    reactant = self.reactant_template.label2atom[label]
                    products_map[a2] = m[reactant]
                else:
                    new_atom = Atom(str(a2))
                    mol_to_transform.atoms.append(new_atom)
                    products_map[a2] = new_atom

            #add all the reaction product bonds:
            for a2 in self.product_template.atoms:
                for neighbor, bond_type in self.product_template.ordered_bonds[a2]:
                    if not neighbor in products_map:
                        # print "missing from mapping!"
                        #create new atom and add to mapping
                        new_atom = Atom(str(neighbor))
                        mol_to_transform.atoms.append(new_atom)
                        products_map[neighbor] = new_atom
                    
                    mol_to_transform.add_bond(products_map[a2], products_map[neighbor], bond_type)
            product_set = (mol_to_transform,)
            out.append(product_set)
        return out

class Retro(object):
    def __init__(self, name, rxn_string_list):
        self.__name__ = name
        self.rxn_list = [Reaction(r) for r in rxn_string_list]
    def __str__(self):
        return self.name
    def run(self, reactants):
        return Retro.run_once(reactants,self.rxn_list)
    @classmethod
    def run_to_completion(cls, reactants, rxn_list):
        products = reactants
        while True:
            reactants = products
            products = self.run_once(reactants, rxn_list)
            #reaction stabilized
            if cls.equal(reactants, products):
                break
        # print "returning..." + str([Chem.MolToSmiles(p) for p in products])
        return products
    @classmethod
    def run_rxn_list(cls, rxn_list, reactants):
        """reactants is a tuple"""
        if not isinstance(reactants,tuple):
            raise
        products = []
        for rxn in rxn_list:
            products += rxn.run(reactants)
        if products:
            return products
        else:
            return (reactants,)
    @classmethod
    def run_once(cls, reactants, rxn_list):
        try:
            products = [p for product_set in cls.run_rxn_list(rxn_list, tuple(reactants)) for p in product_set]
        except ValueError, e:
            #invalid number of reactants
            return reactants
        # for m in products:
        #     print Chem.MolToSmiles(m)
        #     Chem.SanitizeMol(m)
        products = cls.drop_duplicates(products)
        if products:
            return products
        else:
            return reactants
    @classmethod
    def drop_duplicates(cls, mol_list):
        prints = []
        out = []
        for m in mol_list:
            rdkit_mol = Chem.MolFromSmiles(mol_to_smiles(m))
            new_print = FingerprintMols.FingerprintMol(rdkit_mol)
            if any([new_print == p for p in prints]):
                continue
            prints.append(new_print)
            out.append(m)
        return out
    @classmethod
    def equal(cls, mol_list1, mol_list2):
        return len(cls.drop_duplicates(mol_list1+ mol_list2)) == len(mol_list1)

class AtomType(object):
    def __init__(self, name, struct):
        self.name = name
        self.struct = struct
    def bond_summary(self):
        return self.struct.bond_summary(self)
    @property
    def elt(self):
        for e in ["Cl","C","O","Br","I"]:
            if self.name.replace("!","").startswith(e):
                return e
        raise Exception("Couldn't parse AtomPattern " + self.name)
    def __str__(self):
        return self.name
        
class Atom(AtomType):
    """Name, Bonds/Neighbors, Orientation/Stereochemistry, Labels"""
    # def bond_summaries(self):
    #     return [(str(bond_type), str(a)) for bond_type, a in self.bonds]
    def max_bond_cnt(self):
        elt2cnt = {"H":1,
                   "C":4,
                   "O":2,
                   "Br":1,
                   "I":1,
                   "Cl":1}
        return elt2cnt[self.name]
    def cur_bond_cnt(self):
        return len(self.struct.bonds[self])
    def get_hydrogen_cnt(self):
        return self.max_bond_cnt() - self.cur_bond_cnt()

class AtomPattern(AtomType):
    """Contains partial information about an atom 
    for example "!C" (not carbon), or "*" (any atom), "OH" (oxygen attached to >=1 hydrogen)
    """
    def parse_hydrogen_cnt(self):
        """Atom patterns may contain a number of hydrogens, for example
        CH3 or OH
        """
        m = re.findall("H(\d*)",self.name)
        if not m:
            return None
        else:
            if m[0] == "":
                return 1
            else:
                return int(m[0])
    def matches_atom(self, atom):
        if not isinstance(atom, Atom):
            raise
        pattern_hydrogen_cnt = self.parse_hydrogen_cnt()
        atom_hydrogen_cnt = atom.get_hydrogen_cnt()
        if pattern_hydrogen_cnt and not pattern_hydrogen_cnt == atom_hydrogen_cnt:
            return False
        if self.name == "*":
            return True
        elif self.elt == atom.elt:
            return True
        elif self.name.startswith("!") and not atom.name == self.name[1:]:
            return True
        return False
    
    
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

def parse_atom(main):
    if ":" in main:
        return main.split(":",1)
    return main, None
        
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

def find_substructure(substruct, mol):
    #atom -> (bond_type, atom_type) -> cnt
    # mol_cnts = dict((a,counter(a.bond_summaries())) for a in mol.atoms)
    assert(isinstance(substruct, Substructure) and isinstance(mol, Molecule))

    def is_bond_subset(atom1, atom2):
        #need a mapping from atom1's bonds to a subset of atom2's bonds
        bonds1 = atom1.bond_summary()
        n1 = len(bonds1)
        bonds2 = atom2.bond_summary()
        n2 = len(bonds2)
        if n2 < n1:
            return False
        #try all mappings one at a time and see if any work
        for bonds2_subset in itertools.permutations(bonds2,n1):
            match = True
            for b1, b2 in zip(bonds1, bonds2_subset):
                a1, bond_type1 = b1
                a2, bond_type2 = b2
                if not a1.matches_atom(a2) or not bond_type1 == bond_type2:
                    match = False
            if match:
                return True
        return False

    possibilities = {}
    for a1 in substruct.atoms:
        for a2 in mol.atoms:
            if a1.matches_atom(a2) and is_bond_subset(a1, a2):
                possibilities.setdefault(a1,[]).append(a2)

    for a1 in substruct.atoms:
        for a2 in mol.atoms:
            pass

    return match_substructure_helper(substruct, mol, possibilities, {}, {})
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
    m = Molecule(s)
    s_out = mol_to_smiles(m)
    assert(s_out == s)

def test_substructure():
    m1 = Molecule("CCC(Br)CI")
    # m2 = smiles_to_mol('[C:1][C:2]([C:3])(Br)[C:4]') #>>[C:1][C:2]([C:3])=[C:4]')
    m2 = Substructure("C(Br)CI")
    assert(len(find_substructure(m2, m1)) == 1)

    m1 = Molecule("CCC(Br)C")
    # m2 = smiles_to_mol('[C:1][C:2]([C:3])(Br)[C:4]') #>>[C:1][C:2]([C:3])=[C:4]')
    m2 = Substructure("C(Br)C")
    assert(len(find_substructure(m2, m1)) == 2)

def test_reaction():
    m1 = Molecule("CC(C)(Br)CI")
    m2 = Substructure('[C:1][C:2]([C:3])(Br)[C:4]')
    m3 = Substructure('[C:1][C:2]([C:3])=[C:4]')
    rxn = Reaction('[C:1][C:2]([C:3])(Br)[C:4]>>[C:1][C:2]([C:3])=[C:4]')
    print [str(out) for out in rxn.run(m1)]

def test_all():
    test_mol_to_smiles()
    test_substructure()

def test_retro():
    halohydrin_formation = Retro("Halohydrin Formation", ['[C:1][C:2]([C:3])(Br)[C:4][OH]>>[C:1][C:2]([C:3])=[C:4]', '[C:1][C:2](Br)[C:3]([!C:4])([!C:5])[OH]>>[C:1][C:2]=[C:3]([*:4])([*:5])'])
    print halohydrin_formation.run([Molecule("CC(C)(Br)COC")])
    
if __name__ == "__main__":
    # m1 = smiles_to_mol("CCC(Br)CI")
    # print mol_to_smiles(m1)
    test_retro()

    # test_substructure()
