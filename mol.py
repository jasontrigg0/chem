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
        self.setup()
    def setup(self):
        self.atoms = []
        self.atom2label = {} #for labels like [C:4]
        self.label2atom = {}

        #maintain two lists of bonds -- self.bonds for quick access
        #and self.ordered_bonds to store stereochemical information
        self.bonds = {} #atom1 -> atom2 -> bond_type (bonds are symmetric)
        self.ordered_bonds = {} #atom1 -> [(atom2,bond_type)]

        #when the ring_bonds are removed,
        #the rest of the nodes form a tree (forest)
        self.children = {}
        self.parent = {}
        self.ring_bonds = {} #used by DFS: {atom1: [(bond_type,label_number)]}
        self.chirality = {}
        self.atom2component = {} #atom : component#
        self.component2root = {} #component# : atom
        
        self.ring_bond_labels = {} #{label: [(atom,bond_type)]}
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
                bond_string, main, ring_bond_label = re.findall(r"^([{Bond.all_bond_regex}]?)\[(.*)\]([{Bond.all_bond_regex}]?\d*)$".format(**globals()),chunk)[0]
                #add
            else:
                bond_string, main, ring_bond_label = re.findall(r"^([{Bond.all_bond_regex}]?)(.*?)([{Bond.all_bond_regex}]?\d*)$".format(**globals()),chunk)[0]
            print "bond string: " , bond_string
            atom_string, label = parse_atom(main)
            if template:
                next_atom = AtomPattern(atom_string, self)
            else:
                next_atom = Atom(atom_string, self)
            if "@@" in atom_string:
                self.chirality[next_atom] = "@@"
            elif "@" in atom_string:
                self.chirality[next_atom] = "@"
            if label: #C:1
                self.label2atom[label] = next_atom
                self.atom2label[next_atom] = label

            if bond_string in Bond.valid_bond_strings:
                b = Bond(bond_string)
                if atom:
                    self.add_bond(atom, next_atom,b)

            #add ring_bond_label
            if ring_bond_label:
                self.ring_bond_labels.setdefault(ring_bond_label,[]).append((next_atom, b))
                if len(self.ring_bond_labels[ring_bond_label]) > 2:
                    raise Exception("ERROR: three or more copies of the same ring bond")
                elif len(self.ring_bond_labels[ring_bond_label]) == 2:
                    a1,b1 = self.ring_bond_labels[ring_bond_label][0]
                    a2,b2 = self.ring_bond_labels[ring_bond_label][1]
                    if b1 != b2:
                        raise Exception("ERROR: ring bonds with same label ({ring_bond_label}) have different bonds".format(**vars()))
                    
                    self.add_bond(a1,a2,b1)
            self.atoms.append(next_atom)
            atom = self.last_atom()
    def copy(self):
        #return
        #copied_object, along with a map from old atoms to new atoms
        out, atom_map = self.copy_with_map()
        return out
    def copy_with_map(self):
        out = object.__new__(self.__class__)
        out.setup()

        #create a mapping from atoms in the current obj to atoms in the copied object
        atom_map = {}
        for a in self.atoms:
            a_new = Atom(a.name, out)
            atom_map[a] = a_new

        out.atoms = atom_map.values()
            
        out.atom2label = dict((atom_map[a],l) for a,l in self.atom2label.items())
        out.label2atom = dict((l,atom_map[a]) for l,a in self.label2atom.items())

        out.chirality = dict((atom_map[a],c) for a,c in self.chirality.items())
        
        #the self.bonds and self.ordered_bonds dictionaries
        #are nested, so copy.copy won't work
        #also deepcopy will copy the atom elements themselves, which
        #we don't want
        #so re-add the bonds:
        for a,v in self.bonds.items():
            for a2,bond in v.items():
                out.add_bond(atom_map[a],atom_map[a2],bond)


        for a1,l in self.children.items():
            out.children[atom_map[a1]] = [atom_map[a2] for a2 in l]

        out.parent = dict((atom_map[a1],atom_map[a2]) for a1,a2 in self.parent.items())

        for a1,v in self.ring_bonds.items():
            out.ring_bonds.setdefault(atom_map[a1],[]).append(v)
        return out, atom_map
    
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
        if atom1 in self.chirality:
            self.chirality.pop(atom1)
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
        self.atom2component = {}
        self.component2root = {}
        curr_component = 1
        curr_ring_label = 1
        done = set()
        stack = collections.deque()
        while len(self.atom2component) < len(self.atoms):
            #find next atom to start dfs from
            for a in self.atoms:
                if a not in self.atom2component:
                    if self.component2root:
                        curr_component = max(self.component2root.keys()) + 1
                    self.component2root[curr_component] = a
                    self.atom2component[a] = curr_component
                    stack.append((a,None,None))
                    break
            #run dfs from that node
            while stack:
                curr_atom, parent, bond_type = stack.pop()
                if curr_atom in done:
                    #second time we've reached this atom --> ring bond!
                    self.ring_bonds.setdefault(curr_atom,[]).append((bond_type, curr_ring_label))
                    self.ring_bonds.setdefault(parent,[]).append((bond_type, curr_ring_label))
                    curr_ring_label += 1
                else:
                    if parent:
                        self.parent[curr_atom] = parent
                        self.children.setdefault(parent,[]).append(curr_atom)

                    #expand all edges exactly once
                    for neighbor, bond_type in self.ordered_bonds.get(curr_atom,[]):
                        if neighbor in done: continue #don't process the same edge twice
                        stack.append((neighbor, curr_atom, bond_type))
                    self.atom2component[curr_atom] = curr_component
                    done.add(curr_atom)
    def dfs(self):
        for start_node in self.atom2component:
            stack = collections.deque()
            stack.append(start_node)
            while stack:
                next_atom = stack.pop()
                for a in self.children.get(next_atom,[]):
                    stack.append(a)
                yield next_atom
    def __str__(self):
        return mol_to_smiles(self)
    def bond_str(self):
        out = []
        for a1,v in self.bonds.items():
            for a2,b2 in v.items():
                out.append([str(a1),str(a2),str(b2)])
        return str(out)
    def bond_summary_str(self):
        out = ""
        for a1 in self.atoms:
            out += str([(str(x),str(y)) for x,y in self.bond_summary(a1)]) + "\n"
        return out
    
class Molecule(Structure):
    """
    Molecule has all the features of a Structure along with implied Hydrogen atoms
    """
    def __init__(self,smiles):
        super(Molecule,self).__init__()
        #parse the smiles clause and add atom-by-atom
        self.append_smiles_clause(smiles, template=False)
    def __eq__(self, other):
        return (len(self.atoms) == len(other.atoms)) and (find_substructure(self, other) != [])
    def __ne__(self, other):
        return not self == other
    
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
        self.smarts = smarts
        #http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
        #4.5 Component-level grouping of SMARTS
        #C.(C.C)
        smiles1, smiles2 = smarts.split(">>",1)
        smiles1 = smiles1.strip(); smiles2 = smiles2.strip()
        self.reactant_template = Substructure(smiles1)
        self.product_template = Substructure(smiles2)
    def run(self, reactants_tuple):
        reactants_mol = reactants_tuple[0]
        mappings = find_substructure(self.reactant_template, reactants_mol)
        out = []
        for m_basic in mappings:
            mol_to_transform, copy_map = reactants_mol.copy_with_map()

            #m_basic is a map from self.reactant_template -> reactants_mol
            #create m a map from self.reactant_template -> mol_to_transform
            m = dict((k,copy_map[v]) for k,v in m_basic.items())
            
            #remove all the reactant atom bonds:
            for a1 in self.reactant_template.atoms:
                for neighbor, bond_type in self.reactant_template.ordered_bonds[a1]:
                    mol_to_transform.remove_bond(m[a1], m[neighbor], bond_type)
                if not a1 in self.reactant_template.atom2label:
                    mol_to_transform.remove_atom(m[a1])


            def add_atom_copy(old_atom, new_mol):
                old_mol = old_atom.struct
                new_atom = Atom(str(old_atom), new_mol)
                new_mol.atoms.append(new_atom)
                if old_atom in old_mol.chirality:
                    new_mol.chirality[new_atom] = old_mol.chirality[old_atom]
                return new_atom
                
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
                    new_atom = add_atom_copy(a2, mol_to_transform)
                    products_map[a2] = new_atom

            #add all the reaction product bonds:
            for a2 in self.product_template.atoms:
                for neighbor, bond_type in self.product_template.ordered_bonds[a2]:
                    if not neighbor in products_map:
                        # print "missing from mapping!"
                        #create new atom and add to mapping
                        new_atom = add_atom_copy(neighbor, mol_to_transform)
                        products_map[neighbor] = new_atom
                    mol_to_transform.add_bond(products_map[a2], products_map[neighbor], bond_type)
            product_set = (mol_to_transform,)
            out.append(product_set)
        return out
    def __str__(self):
        return self.smarts

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
        return len(cls.drop_duplicates(mol_list1 + mol_list2)) == len(mol_list1)

class AtomType(object):
    def __init__(self, name, struct):
        self.name = name.replace("@","")
        self.struct = struct
    def bond_summary(self):
        return self.struct.bond_summary(self)
    @property
    def elt(self):
        for e in ["Cl","C","O","Br","I","F"]:
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
                   "Cl":1,
                   "F":1,
        }
        return elt2cnt[self.name]
    def cur_bond_cnt(self):
        return len(self.struct.bonds[self])
    def get_hydrogen_cnt(self):
        return self.max_bond_cnt() - self.cur_bond_cnt()
    def matches_atom(self, atom):
        return self.name == atom.name

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
    all_bond_regex = "~=#\."
    valid_bond_strings = "~=#"
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
    def __ne__(self, other):
        return not self == other


def take_next(smiles):
    if smiles[0] == "(":
        #TODO: below regex will fail on nested branches
        branch, rest = re.findall(r"^(\(.*?\))(.*)$",smiles)[0]
        return branch, rest

    re_MAIN = r"^([{Bond.all_bond_regex}]?(Cl|Br|O|I|C|F)\d*)(.*)$".format(**globals())
    m_main = re.findall(re_MAIN,smiles)
    if m_main:
        main, _, rest = m_main[0]
        return main, rest

    re_BRACKET = r"^([{Bond.all_bond_regex}]?\[.*?\]\d*)(.*)$".format(**globals())
    m_bracket = re.findall(re_BRACKET, smiles)
    if m_bracket:
        bracket, rest = m_bracket[0]
        return bracket, rest


def get_chunks(smiles):
    while smiles:
        chunk, smiles = take_next(smiles)
        print chunk, smiles
        if not chunk: raise
        yield chunk

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
    component_smiles = [smiles_helper(mol, r) for _,r in mol.component2root.items()]
    return ".".join(component_smiles)

def compare_stereochemistry(bonded_atoms1, bonded_atoms2):
    #return 1 if two lists of bonded atoms are in the same orientation
    #and 0 otherwise
    if len(bonded_atoms1) != 4 or len(bonded_atoms2) != 4:
        raise
    perm = {}
    for i1,a1 in enumerate(bonded_atoms1):
        for i2,a2 in enumerate(bonded_atoms2):
            if a1 == a2:
                perm[i1] = i2
    if not len(perm) == 4:
        raise Exception("ERROR: couldn't match all 4 atoms")
    return 1 * (permutation_parity(perm) == 1)

def permutation_parity(f):
    #f stores permutation in function form
    #easy to code O(n**2) algo
    #http://math.stackexchange.com/a/1553215
    prod = 1
    for i in f:
        for j in f:
            if j <= i:
                continue
            prod *= (f[i] - f[j]) / float(i - j)
    return prod
        
def opposite_chirality(c):
    if c == "@":
        return "@@"
    elif c == "@@":
        return "@"
    else:
        raise
    
def smiles_helper(mol, atom):
    """Compute smiles of a subsection of the molecule rooted at 'atom'"""
    children = mol.children.get(atom,[])
    parent = mol.parent.get(atom,None)
    out = ""
    if parent:
        out += str(mol.bonds[atom][parent])
    #check how the dfs ordering compares with the stereochemistry:
    stereo_str = ""
    if len(mol.bonds.get(atom,[])) == 4 and atom in mol.chirality:
        mol_neighbors = [n for n,_ in mol.ordered_bonds[atom]]
        smiles_neighbors = []
        if mol.parent.get(atom,None):
            smiles_neighbors += [mol.parent[atom]]
        smiles_neighbors += children
        if compare_stereochemistry(mol_neighbors, smiles_neighbors):
            stereo_str = mol.chirality[atom]
        else:
            stereo_str = opposite_chirality(mol.chirality[atom])
    if stereo_str:
        out += "["+str(atom)+stereo_str+"]"
    else:
        out += str(atom)
    for bond_type, label in mol.ring_bonds.get(atom,[]):
        out += str(bond_type) + str(label)
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
            #check chirality if necessary
            if a1 in sub_mol.chirality and matches[a1] in mol.chirality:
                atoms1 = [matches[n] for n,_ in sub_mol.ordered_bonds[a1]]
                atoms2 = [n for n,_ in mol.ordered_bonds[matches[a1]]]
                order_match = compare_stereochemistry(atoms1,atoms2)
                chirality_match = 1 * (sub_mol.chirality[a1] == mol.chirality[matches[a1]])
                if order_match != chirality_match:
                    return []
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
    m_out = Molecule(s_out)
    assert(m_out == m)

def test_substructure():
    m1 = Molecule("CCC(Br)CI")
    m2 = Substructure("C(Br)CI")
    assert(len(find_substructure(m2, m1)) == 1)

    m1 = Molecule("CCC(Br)C")
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


all_retros = {}
    
def rxn_setup():
    all_retros["halohydrin_formation"] = Retro("Halohydrin Formation", ['[C:1][C:2]([C:3])(Br)[C:4][OH]>>[C:1][C:2]([C:3])=[C:4]', '[C:1][C:2](Br)[C:3]([!C:4])([!C:5])[OH]>>[C:1][C:2]=[C:3]([*:4])([*:5])'])

    all_retros["hydrobromination"] = Retro("Hydrobromination", ['[C:1][C:2]([C:3])(Br)[C:4]>>[C:1][C:2]([C:3])=[C:4]','[C:1][C:2](Br)[C:3]([!C:4])([!C:5])>>[C:1][C:2]=[C:3]([*:4])([*:5])'])

    # dehydrohalogenation = Retro("Dehydrohalogenation", ['[C:1]=[C:2] >> [C:1][C:2][F,Cl,Br,I]'])

    all_retros["alcohol_dehydration"] = Retro("Alcohol Dehydration", ['[C:1]=[C:2] >> [C:1][C:2]O'])

    all_retros["alkene_plus_x2"] = Retro("Alkene Plus X2", ['[*:2][C@@:1](F)([*:3])[C@:4](F)([*:6])[*:5] >> [*:2][C:1]([*:3])=[C:4]([*:5])[*:6]'])
    
    all_retros["halohydrin_formation"] = Retro("Halohydrin Formation", ['[C:1][C:2]([C:3])(Br)[C:4][OH]>>[C:1][C:2]([C:3])=[C:4]', '[C:1][C:2](Br)[C:3]([!C:4])([!C:5])[OH]>>[C:1][C:2]=[C:3]([*:4])([*:5])'])

    all_retros["oxymercuration"] = Retro("Oxymercuration", ['[C:1][C:2]([C:3])([OH])[CH:4]>>[C:1][C:2]([C:3])=[C:4]', '[C:1][C:2]([OH])[CH:3]([!C:4])([!C:5])>>[C:1][C:2]=[C:3]([*:4])([*:5])'])

    all_retros["hydroboration"] = Retro("Hydroboration", ['[C:1][C:2]([C:3])[CH:4]([OH])>>[C:1][C:2]([C:3])=[C:4]', '[C:1][C:2][CH2:3][OH]>>[C:1][C:2]=[C:3]', '[C:1][CH:2]([OH])[CH2:3][C:4]>>[C:1][C:2]=[C:3][C:4]'])

    all_retros["alkene_hydrogenation"] = Retro("Alkene Hydrogenation", ['[CH2:1][CH3:2]>>[CH1:1]=[CH2:2]', '[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]', '[CH1:1][CH3:2]>>[CH1:1]=[CH2:2]', '[CH2:1][CH2:2]>>[CH:1]=[CH:2]', '[CH2:1][CH:2]>>[CH:1]=[CH0:2]', '[CH:1][CH:2]>>[CH0:1]=[CH0:2]'])

    all_retros["alkene_hydroxylation"] = Retro("Alkene Hydroxylation", ['[C:1]([OH])[C:2]([OH])>>[C:1]=[C:2]'])

    all_retros["dichlorocarbene_addition"] = Retro("Dichlorocarbene Addition", ['[C:1]C(Cl)(Cl)[C:2]>>[C:1]=[C:2]'])

    all_retros["simmons_smith"] = Retro("Simmons-Smith", ['[C:1]1[CH2][C:2]1>>[C:1]=[C:2]'])

    all_retros["ozonolysis"] = Retro("Ozonolysis", ["[C:1]=O.[C:2]=O >> [C:1]=[C:2]"])
    
def test_retro():
    print "Creating reactant"
    start1 = Molecule("CC(C)(Br)CO")
    print "running retro:"
    a = all_retros["halohydrin_formation"].run([start1])[0]
    print "creating a CC(C)=C molecule"
    b = Molecule("CC(C)=C")
    print "checking substructure between reaction product and copy"
    find_substructure(a,b)

    print "reaction product bonds"
    print a.bond_str()

    print a.bond_summary_str()
    
    print all_retros["halohydrin_formation"].run([start1])[0] == Molecule("CC(C)=C")
    assert(all_retros["halohydrin_formation"].run([start1]) == [Molecule("CC(C)=C")])

    start2 = Molecule("CC(C)(Br)COC")
    assert(all_retros["halohydrin_formation"].run([start2])[0] == start2)

def test_reactions():
    rxn_setup()
    
    #hydrobromination
    start1 = Molecule("CC(C)(Br)C")
    assert(all_retros["hydrobromination"].run([start1])[0] == Molecule("CC(C)=C"))

    #dehydrohalogenation
    # start1 = Molecule("CCBr")

    #alcohol_dehydration
    start1 = Molecule("C=CCCC")
    #TODO:
    #smarter SMILES printing: eg
    #C(CC(C)O)C
    #would look better as
    #CCCC(O)C
    end = [Molecule("CCCCCO"),Molecule("CCCC(O)C")]
    assert(all([output in end for output in all_retros["alcohol_dehydration"].run([start1])]))

    #[*:2][C@@:1](F)([*:3])[C@:4](F)([*:6])[*:5] >> [*:2][C:1]([*:3])=[C:4]([*:5])[*:6]
    start1 = Molecule("C[C@@](F)(Cl)[C@](F)(C)C")
    assert(all_retros["alkene_plus_x2"].run([start1])[0] == Molecule("ClC(C)=C(C)C"))

    start1 = Molecule("CC(C)(Br)CO")
    end1 = Molecule("CC(C)=C")
    assert(all_retros["halohydrin_formation"].run([start1])[0] == end1)

    start1 = Molecule("CC(O)C(Br)(I)")
    end1 = Molecule("C(C)=C(Br)I")
    assert(all_retros["oxymercuration"].run([start1])[0] == end1)

    start1 = Molecule("C(Br)(Br)C")
    end1 = Molecule("C(Br)(Br)=C")
    assert(all_retros["alkene_hydrogenation"].run([start1])[0] == end1)

    start1 = Molecule("BrC(Br)(Br)C(Br)(Br)")
    end1 = start1
    assert(all_retros["alkene_hydroxylation"].run([start1])[0] == end1)
    
    start1 = Molecule("OCCO")
    end1 = Molecule("C=C")
    assert(all_retros["alkene_hydroxylation"].run([start1])[0] == end1)

    start1 = Molecule("CC(Cl)(Cl)C")
    end1 = Molecule("C=C")
    assert(all_retros["dichlorocarbene_addition"].run([start1])[0] == end1)

    start1 = Molecule("CCC1CC1")
    end1 = Molecule("CCC=C")
    assert(all_retros["simmons_smith"].run([start1])[0] == end1)


    
def test_chirality():
    a1 = Molecule("C[C@@](Br)(Cl)I")
    a2 = Molecule("C[C@@](Cl)(Br)I")
    assert(not find_substructure(a1,a2))

    a1 = Molecule("C[C@@](Br)(Cl)I")
    a2 = Molecule("C[C@](Cl)(Br)I")
    assert(find_substructure(a1,a2))

    test_rxn = Reaction("C[C@](Cl)(Br)I >> C[C@@](Cl)(Br)I")
    test_rxn2 = Reaction("C[C@](Br)(Cl)I >> C[C@](Cl)(Br)I")
    r = Molecule("C[C@](Cl)(Br)I")
    out = test_rxn.run((r,))[0][0]
    out2 = test_rxn2.run((out,))[0][0]
    assert(r != out)
    assert(r == out2)
    
def test_smiles_ring_bonds():
    a1 = Molecule("C1CCC1")
    print a1.bond_summary_str()

def test_copy():
    m = Molecule("[C:1]1CCC1")
    m.setup_dfs()
    m2 = m.copy()
    print m2 == m
    
if __name__ == "__main__":
    # print mol_to_smiles(m1)

    # test_reactions()
    # test_smiles_ring_bonds()
    # print list(get_chunks("[C:1]1CCC1"))
    # test_reactions()


    # reactants = []
    # esterification = Reaction("C(=O)O.OCC>>C(=O)OCC.O")
    # print esterification.run(reactants)

    # m= Molecule("C(=O)O.OCC") #>>C(=O)OCC.O")
    # m.setup_dfs()
    # print m

    print Molecule("C.C.C")
