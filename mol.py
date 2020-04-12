#!/usr/bin/env python
import re

import utils
import Queue
import collections
import copy
import itertools

from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem
from rdkit import DataStructs

logger = utils.Logger()


#SMILES spec:
#http://opensmiles.org/spec/open-smiles-3-input.html
#(lots left to implement!)

#NOTE:
#not following standard SMILES for cis-trans stereochemistry because of the issue here:
#(1)http://opensmiles.org/spec/open-smiles-6-extensions.html#6.3
#(2)http://baoilleach.blogspot.com/2010/09/are-you-on-my-side-or-not-its-ez-part.html
#(3)http://forums.openbabel.org/Cis-trans-stereochem-with-SMILES-at-ring-closures-td2715247.html

#Instead using the syntax proposed here:
#http://opensmiles.org/spec/open-smiles-6-extensions.html#6.3

class Structure(object):
    """Representation of a molecular structure:
    """
    def __init__(self):
        self.setup()
    def setup(self):
        self.explicit_atoms = []
        self.implicit_atoms = []

        self.atom2label = {} #for labels like [C:4]
        self.label2atom = {}

        #maintain two lists of bonds -- self.bonds for quick access
        #and self.ordered_bonds to store stereochemical information
        self.bonds = {} #atom1 -> atom2 -> bond_type (bonds are symmetric)
        self.implicit_bonds = {} #implicit hydrogen bonds atom1 -> atom2 -> bond_type

        self.ordered_bonds = {} #atom1 -> [(atom2,bond_type)]

        #component_bonds maintains a list of componant constraints
        #in the case of templates
        #ie atom1, atom2 must be in the same component (1)
        #or atom1, atom2 must be in different components (0)
        self.component_bonds = {} #atom1 -> [(atom2, is_same_component_boolean)]

        #when the ring_bonds are removed,
        #the rest of the nodes form a tree (forest)
        self.children = {}
        self.parent = {}
        self.ring_bonds = {} #used by DFS: {atom1: [(bond_type,label_number)]}
        self.ring_bond_labels = {} #{label: [(atom,bond_type)]}

        #chirality of tetrahedral bonds
        self.chirality = {}

        #store information about the cis-trans nature of C=C bonds
        self.c_double_bonds = [] #[set([atom1,atom2])]
        self.c_double_bond_atoms = set() #any atoms in a C=C bond

        self.atom2component = {} #atom : component#
        self.component2root = {} #component# : atom

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
            if bond_string in Bond.valid_bond_strings: #includes bond_string = ""
                b = Bond(bond_string)
                if atom:
                    self.add_bond(atom, next_atom,b)
            elif re.findall("["+Bond.all_bond_regex+"]",bond_string): #bond is "." or "`" or "'"
                if atom:
                    self.add_component_bond(atom, next_atom, bond_string)

            #add ring_bond_label
            if ring_bond_label:
                self.ring_bond_labels.setdefault(ring_bond_label,[]).append((next_atom, b))
                if len(self.ring_bond_labels[ring_bond_label]) > 2:
                    raise Exception("ERROR: three or more copies of the same ring bond label")
                elif len(self.ring_bond_labels[ring_bond_label]) == 2:
                    a1,b1 = self.ring_bond_labels[ring_bond_label][0]
                    a2,b2 = self.ring_bond_labels[ring_bond_label][1]
                    if b1 != b2:
                        raise Exception("ERROR: ring bonds with same label ({ring_bond_label}) have different bonds".format(**vars()))
                    self.add_bond(a1,a2,b1)
            self.add_atom(next_atom)
            atom = self.last_atom()
    def add_implicit_hydrogens(self):
        for a in self.explicit_atoms:
            if a.implicit_hydrogen_cnt():
                for i in range(a.implicit_hydrogen_cnt()):
                    h_atom = Atom("H",self)
                    self.implicit_atoms.append(h_atom)
                    self.implicit_bonds.setdefault(h_atom,{})[a] = Bond("")
                    self.implicit_bonds.setdefault(a,{})[h_atom] = Bond("")
                    #TODO check that we're doing everything in self.add_bond
                    # self.add_bond(a,h_atom,Bond(""))
    @classmethod
    def merge(cls, list_of_structures):
        """
        combine multiple structures into
        one structure with multiple components
        """
        out = object.__new__(cls)
        out.setup()
        out.explicit_atoms = [a for s in list_of_structures for a in s.explicit_atoms]
        out.implicit_atoms = [a for s in list_of_structures for a in s.implicit_atoms]
        for s in list_of_structures:
            out.chirality.update(s.chirality)
            if s.atom2label:
                raise
            if s.label2atom:
                raise
            out.bonds.update(s.bonds)
            out.ordered_bonds.update(s.ordered_bonds)
            out.implicit_bonds.update(s.implicit_bonds)

            out.c_double_bonds += s.c_double_bonds
            out.c_double_bond_atoms.update(s.c_double_bond_atoms)
        
        out.setup_dfs()
        return out
    def split(self):
        """
        Return list of components as individual structures
        """
        self.setup_dfs()
        out = []
        for component in self.component2root:
            obj = object.__new__(self.__class__)
            obj.setup()
            obj.explicit_atoms = [a for a in self.explicit_atoms if self.atom2component[a] == component]
            component_atoms = set(obj.explicit_atoms)
            obj.atom2label = dict((a,l) for a,l in self.atom2label.items() if a in component_atoms)
            obj.label2atom = dict((l,a) for l,a in self.label2atom.items() if a in component_atoms)
            obj.chirality = dict((a,c) for a,c in self.chirality.items() if a in component_atoms)
            obj.bonds = dict((a,v) for a,v in self.bonds.items() if a in component_atoms)
            obj.ordered_bonds = dict((a,v) for a,v in self.ordered_bonds.items() if a in component_atoms)
            obj.c_double_bonds = [set([a1,a2]) for a1,a2 in self.c_double_bonds if a1 in component_atoms and a2 in component_atoms]
            obj.c_double_bond_atoms = set([a for a in self.c_double_bond_atoms if a in component_atoms])
        
            for a in obj.explicit_atoms:
                a.struct = obj
                
            obj.setup_dfs()
            if isinstance(obj, Molecule):
                obj.add_implicit_hydrogens()
            out.append(obj)
        return out
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
        for a in self.explicit_atoms + self.implicit_atoms:
            a_new = Atom(a.name, out)
            atom_map[a] = a_new
            
        out.explicit_atoms = [atom_map[a] for a in self.explicit_atoms]
        out.implicit_atoms = [atom_map[a] for a in self.implicit_atoms]
            
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
                try:
                    out.add_bond(atom_map[a],atom_map[a2],bond)
                except:
                    deep_print(self.explicit_atoms)
                    deep_print(self.implicit_atoms)
                    deep_print(self.bonds)
                    print a,a2
                    print self
                    raise

        out.c_double_bonds = [set([atom_map[a1],atom_map[a2]]) for a1,a2 in self.c_double_bonds]
        out.c_double_bond_atoms = set([atom_map[a] for a in self.c_double_bond_atoms])

        out.setup_dfs() #sets up the ring bonds
        return out, atom_map
    
    def last_atom(self):
        if self.explicit_atoms:
            return self.explicit_atoms[-1]
        else:
            return None
    def add_atom(self, atom):
        self.explicit_atoms.append(atom)
    def add_bond(self, atom1, atom2, bond_type):
        if self.bonds.get(atom1,{}).get(atom2,None):
            return
        self.bonds.setdefault(atom1,{})[atom2] = bond_type
        self.bonds.setdefault(atom2,{})[atom1] = bond_type

        self.ordered_bonds.setdefault(atom1,[]).append((atom2,bond_type))
        self.ordered_bonds.setdefault(atom2,[]).append((atom1,bond_type))

        if bond_type == Bond("=") and str(atom1) == "C" and str(atom2) == "C":
            self.c_double_bonds.append(set([atom1,atom2]))
            self.c_double_bond_atoms.add(atom1)
            self.c_double_bond_atoms.add(atom2)
                        
    def add_component_bond(self, atom1, atom2, bond_string):
        if bond_string == INTER_CHAR:
            is_same_component_boolean = 0
        elif bond_string == INTRA_CHAR:
            is_same_component_boolean = 1
        else:
            return
        self.component_bonds.setdefault(atom1,[]).append((atom2, is_same_component_boolean))
    def add_ring_bond(self, atom1, atom2, bond_type):
        if self.ring_bond_labels:
            next_ring_label = max([int(i) for i in self.ring_bond_labels.keys()]) + 1
        else:
            next_ring_label = 1
        self.ring_bonds.setdefault(atom1,[]).append((bond_type, next_ring_label))
        self.ring_bonds.setdefault(atom2,[]).append((bond_type, next_ring_label))
        self.ring_bond_labels.setdefault(next_ring_label,[]).append((atom1, bond_type))
        self.ring_bond_labels.setdefault(next_ring_label,[]).append((atom2, bond_type))
    def get_bonds(self, atom):
        out = []
        for a2, bond_type in self.bonds.get(atom,{}).items() + self.implicit_bonds.get(atom,{}).items():
            out.append((a2, bond_type))
        return out
    # def bond_summary(self, atom):
    #     out = []
    #     for a2, bond_type in self.ordered_bonds.get(atom,[]):
    #         out.append((a2, bond_type))
    # return out
    def all_bonds(self, raw=False):
        out = []
        for a1 in self.explicit_atoms:
            for a2, bond_type in self.ordered_bonds.get(a1,[]):
                if raw:
                    out.append((a1, a2, bond_type))
                else:
                    out.append((str(a1), str(a2), str(bond_type)))
        return out
    def remove_atom(self, atom1):
        if atom1 in self.explicit_atoms:
            self.explicit_atoms.remove(atom1)
            self.bonds.pop(atom1)
            self.ordered_bonds.pop(atom1)
        else:
            self.implicit_atoms.remove(atom1)
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
        if bond_type == Bond("=") and str(atom1) == "C" and str(atom2) == "C":
            self.c_double_bonds.remove(set([atom1,atom2]))
            
    def setup_dfs(self):
        self.children = {}
        self.parent = {}
        self.ring_bonds = {}
        self.ring_bond_labels = {}
        self.atom2component = {}
        self.component2root = {}
        curr_component = 1
        done = set()
        stack = collections.deque()
        while len(self.atom2component) < len(self.explicit_atoms):
            #find next atom to start dfs from
            #start with one of the atoms with the lowest degree to make prettier smiles represenation
            unused_atoms = sorted([a for a in self.explicit_atoms if a not in self.atom2component],key=lambda x: len([a for a in self.bonds.get(x,{}) if a.elt != "H"])) 
            next_atom = unused_atoms[0]
            if self.component2root:
                curr_component = max(self.component2root.keys()) + 1
            self.component2root[curr_component] = next_atom
            self.atom2component[next_atom] = curr_component
            stack.append((next_atom,None,None))
            #run dfs from that node
            while stack:
                curr_atom, parent, bond_type = stack.pop()
                if curr_atom in done:
                    #second time we've reached this atom --> ring bond!
                    self.add_ring_bond(curr_atom, parent, bond_type)
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
        for a1 in self.explicit_atoms:
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
        self.add_implicit_hydrogens()
    def __eq__(self, other):
        if not type(self) == type(other):
            return False
        return (len(self.explicit_atoms) == len(other.explicit_atoms)) and (find_substructure(self, other) != [])
    def __ne__(self, other):
        return not self == other
    def __hash__(self):
        rdkit_mol = Chem.MolFromSmiles(str(self))
        new_print = FingerprintMols.FingerprintMol(rdkit_mol)
        num_bits = new_print.GetNumBits()
        bits = "".join([str(1 *new_print.GetBit(i)) for i in range(num_bits)])
        return int(bits,2)
    def similarity(self, other):
        return self.rdkit_sim(other) #0.5 * (self.rdkit_sim(other) + self.autochem_sim(other))
    def rdkit_sim(self, other):
        if not isinstance(other, Molecule):
            raise
        rdkit1 = Chem.MolFromSmiles(str(self))
        f_end = FingerprintMols.FingerprintMol(rdkit1)
        rdkit2 = Chem.MolFromSmiles(str(other))
        f_start = FingerprintMols.FingerprintMol(rdkit2)
        sim = DataStructs.FingerprintSimilarity(f_end, f_start)
        return sim
    def autochem_sim(self, other):
        """create our own measure of similarity because the rdkit measure 
        isn't very effective for searching
        """
        #carbon cnt
        c1 = len([a for a in self.explicit_atoms if a.elt == "C"])
        c2 = len([a for a in other.explicit_atoms if a.elt == "C"])
        carbon_sim = 1 * (c1 == c2)
        
        #bonds on the carbon backbone
        def get_cbond_cnts(m):
            cbonds1 = []
            for a1 in m.bonds:
                if not a1.elt == "C": continue
                for a2 in m.bonds[a1]:
                    if not a2.elt == "C": continue
                    cbonds1.append(m.bonds[a1][a2])
            cbonds_cnt = counter(cbonds1, key = lambda x: x.bond_cnt)
            return cbonds_cnt

        cbonds1_cnt = get_cbond_cnts(self)
        cbonds2_cnt = get_cbond_cnts(other)
        matched = 0
        total = 0
        for bond_type in set(cbonds1_cnt.keys() + cbonds2_cnt.keys()):
            cnt1 = cbonds1_cnt.get(bond_type,0)
            cnt2 = cbonds2_cnt.get(bond_type,0)
            matched += min(cnt1, cnt2)
            total += max(cnt1, cnt2)
        single = max(cbonds1_cnt.get(1,0), cbonds2_cnt.get(2,0))
        non_single = 1 - single / float(total) #pct of carbon bonds that are not single bonds
        cbond_sim = matched / float(total)

        #penalize complexity somewhat
        total = 0.5 * (carbon_sim + cbond_sim)
        # total -= (c1 != c2) * (2 - 1 / float(c1) - 1 / float(c2))
        # total -= (cbond_sim != 1) * (non_single)
        
        #TODO: non-carbon cnts
        return total
        
class Substructure(Structure):
    """
    """
    def __init__(self,smiles,constraints=None):
        super(Substructure,self).__init__()
        smiles = self.preprocess_smiles(smiles)
        self.append_smiles_clause(smiles, template=True)

        self.constraints = []
        if constraints:
            self.constraints = constraints
    @classmethod
    def preprocess_smiles(self, smiles):
        #substructures have special meaning for "."
        #in different contexts
        #preprocess the smiles string for these cases

        #list of which locations represent specifically bonds within or between components
        inter = []
        intra = []
        for i,c in enumerate(smiles):
            if c == ".":
                left_cnt = len([c2 for c2 in smiles[:i] if c2 == "("])
                right_cnt = len([c2 for c2 in smiles[:i] if c2 == ")"])
                if left_cnt > right_cnt:
                    intra.append(i)
                elif smiles[i-1] == ")" and smiles[i+1] == "(":
                    inter.append(i)

        #remove the top-level parentheses
        depth = 0
        to_delete = []
        for i,c in enumerate(smiles):
            if depth == 0 and c == "(":
                if i > 0 and smiles[i-1] == ".":
                    to_delete.append(i)
            if depth == 1 and c == ")":
                if i < (len(smiles)-1) and smiles[i+1] == ".":
                    to_delete.append(i)
            if c == "(":
                depth += 1
            elif c == ")":
                depth -= 1

        for i in inter:
            smiles = smiles[:i] + INTER_CHAR + smiles[i+1:]
        for i in intra:
            smiles = smiles[:i] + INTRA_CHAR + smiles[i+1:]
        smiles = [s for i,s in enumerate(smiles) if i not in to_delete]
        smiles = "".join(smiles)
        return smiles

class Constraint(object):
    def test(self, struct, mol, matches):
        raise Exception("Unimplemented")

class Markovnikov(Constraint):
    def __init__(self, label1, label2):
        """initialize markovnikov with -- the carbon at label1 is *more*
        substituted than the carbon at label2
        """
        self.label1 = str(label1)
        self.label2 = str(label2)
    def test(self, struct, mol, matches):
        if not self.label1 in struct.label2atom or not self.label2 in struct.label2atom:
            raise Exception("ERROR: markovnikov labels not present in substructure")
        a1 = matches[struct.label2atom[self.label1]]
        a2 = matches[struct.label2atom[self.label2]]
        if not a1.elt == "C" or not a2.elt == "C":
            raise Exception("ERROR: non-carbon markovnikov constraint atoms detected")
        
        subs1 = len([n for n,_ in mol.ordered_bonds[a1] if n.elt == "C"])
        subs2 = len([n for n,_ in mol.ordered_bonds[a2] if n.elt == "C"])
        return subs1 > subs2
        
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
    def __init__(self, smarts, constraints = None):
        self.smarts = smarts
        #http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
        #4.5 Component-level grouping of SMARTS
        #C.(C.C)
        smiles1, smiles2 = smarts.split(">>",1)
        smiles1 = smiles1.strip(); smiles2 = smiles2.strip()
        self.reactant_template = Substructure(smiles1, constraints = constraints)
        self.product_template = Substructure(smiles2)
    def run(self, reactants_tuple):
        reactants_mol = Molecule.merge(reactants_tuple)
        mappings = find_substructure(self.reactant_template, reactants_mol)
        out = []
        for m_basic in mappings:
            mol_to_transform, copy_map = reactants_mol.copy_with_map()
            #m_basic is a map from self.reactant_template -> reactants_mol
            #create m a map from self.reactant_template -> mol_to_transform
            m = dict((k,copy_map[v]) for k,v in m_basic.items())
            #remove all the reactant atom bonds:
            for a1 in self.reactant_template.explicit_atoms:
                for neighbor, bond_type in self.reactant_template.ordered_bonds[a1]:
                    mol_to_transform.remove_bond(m[a1], m[neighbor], bond_type)
                if not a1 in self.reactant_template.atom2label:
                    mol_to_transform.remove_atom(m[a1])

            def add_atom_copy(old_atom, new_mol):
                old_mol = old_atom.struct
                new_atom = Atom(str(old_atom), new_mol)
                new_mol.add_atom(new_atom)
                if old_atom in old_mol.chirality:
                    new_mol.chirality[new_atom] = old_mol.chirality[old_atom]
                return new_atom
                
            #generate mapping from product_template to reactants_mol using
            #the labels to map from product_template to reactant_template + 
            #the existing mapping from reactant_template to reactants_mol
            products_map = {}
            for a2 in self.product_template.explicit_atoms:
                if a2 in self.product_template.atom2label:
                    label = self.product_template.atom2label[a2]
                    reactant = self.reactant_template.label2atom[label]
                    products_map[a2] = m[reactant]
                    # below line prints the label mappings (eg Cl -> 1, Br -> 2, ...)
                    # print label, products_map[a2]
                else:
                    new_atom = add_atom_copy(a2, mol_to_transform)
                    products_map[a2] = new_atom

            #add all the product template bonds:
            for a2 in self.product_template.explicit_atoms:
                for neighbor, bond_type in self.product_template.ordered_bonds.get(a2,[]):
                    if not neighbor in products_map:
                        new_atom = add_atom_copy(neighbor, mol_to_transform)
                        products_map[neighbor] = new_atom
                    if products_map[a2] in mol_to_transform.explicit_atoms and products_map[neighbor] in mol_to_transform.explicit_atoms: #bugfix: only add bonds between explicit atoms in reactant
                        mol_to_transform.add_bond(products_map[a2], products_map[neighbor], bond_type)
                        
            #update mol_to_transform.chirality:
            for a in self.product_template.chirality:
                mol_to_transform.chirality[products_map[a]] = self.product_template.chirality[a]
            for a in mol_to_transform.explicit_atoms:
                if a in mol_to_transform.chirality and not a.is_chiral():
                    mol_to_transform.chirality.pop(a)

            #remove all implicit atoms+bonds:
            mol_to_transform.implicit_atoms = []
            mol_to_transform.implicit_bonds = {}
            #add back implicit atoms + bonds
            mol_to_transform.add_implicit_hydrogens()
            product_set = tuple(mol_to_transform.split(),)
            out.append(product_set)
        return out
    def __str__(self):
        return self.smarts

class Retro(object):
    def __init__(self, name, rxn_list):
        self.__name__ = name
        self.rxn_list = rxn_list
    def __str__(self):
        return self.__name__
    def run(self, reactants):
        if not (isinstance(reactants,tuple) or isinstance(reactants,list)):
            raise Exception("Reactants not a tuple/list type")
        if not all([isinstance(r,Structure) for r in reactants]):
            raise
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
        product_sets = cls.run_rxn_list(rxn_list, tuple(reactants))
        product_sets = list(set(product_sets))
        if product_sets:
            return product_sets
        else:
            return [tuple(reactants)]

class AtomType(object):
    def __init__(self, name, struct):
        self.name = name.replace("@","")
        self.struct = struct
    def bond_summary(self):
        return self.struct.get_bonds(self)
    @property
    def elt(self):
        for e in ["Cl","C","O","Br","I","F","H"]:
            if self.name.replace("!","").startswith(e): #TODO: why replace "!" with "" ?????
                return e
        return None
    def implicit_hydrogen_cnt(self):
        return self.max_bond_cnt() - self.cur_bond_cnt()
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
        return sum(b.bond_cnt for b in self.bonds().values())
    def cur_neighbor_cnt(self):
        return len(self.bonds())
    def is_chiral(self):
        bond_cnt = self.cur_bond_cnt()
        neighbor_cnt = self.cur_neighbor_cnt()
        if neighbor_cnt <= 2:
            return False #carbons not chiral in BrC#CBr
        return bond_cnt == 4 or (bond_cnt == 3 and self.elt == "C")
    def bonds(self):
        return self.struct.bonds.get(self,{})
    def __str__(self):
        return self.name
        
class Atom(AtomType):
    """Name, Bonds/Neighbors, Orientation/Stereochemistry, Labels"""
    def __init__(self, name, struct):
        super(Atom,self).__init__(name, struct)
        self.name = re.sub(r"(.+)H\d*",r"\1",self.name) #OH -> O
    def hydrogen_bond_cnt(self):
        return len([a for a in self.bonds().keys() if a.elt == "H"])
    def hydrogen_cnt(self):
        return self.max_bond_cnt() - self.cur_bond_cnt() + self.hydrogen_bond_cnt()
    def matches_atom(self, atom):
        return self.name == atom.name

class AtomPattern(AtomType):
    """Contains partial information about an atom 
    for example "!C" (not carbon), or "*" (any atom), "OH" (oxygen attached to >=1 hydrogen)
    """
    def hydrogen_cnt(self):
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
        if self.name == "*":
            return True
        pattern_hydrogen_cnt = self.hydrogen_cnt()
        atom_hydrogen_cnt = atom.hydrogen_cnt()
        if pattern_hydrogen_cnt  and not pattern_hydrogen_cnt == atom_hydrogen_cnt:
            return False
        elif self.elt == atom.elt:
            return True
        elif self.name.startswith("!") and not atom.name == self.name[1:]:
            return True
        return False
    

INTER_CHAR = "`"
INTRA_CHAR = "'"
    
class Bond(object):
    valid_bond_strings = r"~=#"
    all_bond_regex = valid_bond_strings + "\." + INTER_CHAR + INTRA_CHAR
    def __init__(self, bond_string):
        if bond_string in ("\\","/"):
            bond_string = ""
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
        depth=0
        for i,x in enumerate(smiles):
            if x == "(": depth += 1
            if x == ")": depth -= 1
            if depth == 0:
                return smiles[:(i+1)], smiles[(i+1):]
        raise
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

    raise Exception("Unable to parse smiles substring: {smiles}".format(**vars()))

def get_chunks(smiles):
    while smiles:
        chunk, smiles = take_next(smiles)
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
    if not mol.explicit_atoms:
        return ""
    component_smiles = [smiles_helper(mol, r) for _,r in mol.component2root.items()]
    return ".".join(component_smiles)

def compare_stereochemistry(bonded_atoms1, bonded_atoms2):
    #return 1 if two lists of bonded atoms are in the same orientation
    #and 0 otherwise
    if not len(bonded_atoms1) == len(bonded_atoms2):
        raise
    allenal = (len(bonded_atoms1) == 3 and len(bonded_atoms2) == 3)
    tetrahedral = (len(bonded_atoms1) == 4 and len(bonded_atoms2) == 4)
    if not allenal and not tetrahedral:
        deep_print(bonded_atoms1)
        deep_print(bonded_atoms2)
        raise Exception("trying to compare stereochemistry of two atoms with unequal neighbor count or < 3 neighbors (ie non-chiral)")
    perm = {}
    for i1,a1 in enumerate(bonded_atoms1):
        for i2,a2 in enumerate(bonded_atoms2):
            if a1 == a2:
                perm[i1] = i2
    if len(perm) != len(bonded_atoms1):
        raise Exception("Couldn't match all atoms")
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
    children = [a for a in mol.children.get(atom,[]) if a.elt != "H"]
    children = sorted(children,key=lambda x: len(mol.bonds.get(x,{})))
    parent = mol.parent.get(atom,None)
    out = ""
    if parent:
        out += str(mol.bonds[atom][parent])
    #check how the dfs ordering compares with the stereochemistry:
    stereo_str = ""
    if atom in mol.chirality and atom.is_chiral():
        mol_neighbors = [n for n,_ in mol.ordered_bonds[atom]]
        smiles_neighbors = []
        if parent:
            smiles_neighbors += [parent]
        if atom in mol.ring_bonds:
            #Currently only supporting one ring bond for an atom --
            #is >1 possible?
            assert(len(mol.ring_bonds[atom]) == 1)
            label = mol.ring_bonds[atom][0][1]
            assert(len(mol.ring_bond_labels[label]) == 2)
            smiles_neighbors += [i[0] for i in mol.ring_bond_labels[label] if i[0] != atom]
        smiles_neighbors += children
        if atom.implicit_hydrogen_cnt() > 0:
            assert(len(smiles_neighbors) == len(mol_neighbors))
            #implicit hydrogen is always in the second position
            h_atom = Atom("H", None)
            mol_neighbors = mol_neighbors[:1] + [h_atom] + mol_neighbors[1:]
            smiles_neighbors = smiles_neighbors[:1] + [h_atom] + smiles_neighbors[1:]
        if compare_stereochemistry(mol_neighbors, smiles_neighbors):
            stereo_str = mol.chirality[atom]
        else:
            stereo_str = opposite_chirality(mol.chirality[atom])
    atom_str = str(atom)
    if atom_str != "H":
        if stereo_str or (atom_str != atom.elt):
            out += "["+atom_str+stereo_str+"]"
        else:
            out += atom_str
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

# def is_subset(list1, list2):
#     return is_counter_subset(counter(list1),counter(list2))

# def is_counter_subset(counter1, counter2):
#     #check that counter1 < counter2
#     for k1,v1 in counter1.items():
#         if v1 > counter2.get(k1,0):
#             return False
#     return True


########substructure matching########
def find_substructure(substruct, mol):
    #atom -> (bond_type, atom_type) -> cnt
    # mol_cnts = dict((a,counter(a.bond_summaries())) for a in mol.atoms)
    substruct.setup_dfs()
    mol.setup_dfs()

    possibilities = gen_possibilities(substruct, mol, {})
    return match_substructure_helper(substruct, mol, possibilities, {}, {})

def is_bond_subset(atom1, atom2, matches):
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
            if not a1.matches_atom(a2) or not bond_type1 == bond_type2 or (a1 in matches and a2 != matches[a1]):
                match = False
        if match:
            # print "true"
            return True
    # print "false"
    return False

def gen_possibilities(sub_mol, mol, matches):
    possibilities = {}
    for a1 in sub_mol.explicit_atoms:
        for a2 in mol.explicit_atoms + mol.implicit_atoms: #allow matching to hydrogens if necessary
            if a1.matches_atom(a2) and is_bond_subset(a1, a2, matches):
                possibilities.setdefault(a1,[]).append(a2)
    return possibilities

def match_substructure_helper(sub_mol, mol, possibilities, matches, done):
    """matches = [(atom from sub_mol, mapped atom from mol)]"""
    #base case
    #all done matching -- check that it works!
    if len(matches) == len(sub_mol.explicit_atoms):
        for a1 in sub_mol.explicit_atoms:
            b1 = matches[a1]
            if b1 == None:
                continue
            for neighbor, bond_type in sub_mol.ordered_bonds.get(a1,[]):
                if neighbor.elt == "H": continue #only worry about backbone bonds
                if matches[neighbor] is None: continue
                if not (matches[neighbor], bond_type) in mol.get_bonds(b1):
                    return [] #no match!
            #check constraints (eg markovnikov)
            if isinstance(sub_mol,Substructure):
                for c in sub_mol.constraints:
                    if not c.test(sub_mol, mol, matches):
                        return []
        for a1 in sub_mol.explicit_atoms:
            #check tetrahedral chirality if necessary:
            if a1 in sub_mol.chirality and matches[a1] in mol.chirality and a1 not in sub_mol.c_double_bond_atoms and matches[a1] not in mol.c_double_bond_atoms:
                atoms1 = [matches[n] for n,_ in sub_mol.ordered_bonds[a1]]
                atoms2 = [n for n,_ in mol.ordered_bonds[matches[a1]]]
                #treat implicit hydrogens as the second bonded atom:
                #http://opensmiles.org/spec/open-smiles-3-input.html
                # deep_print(sub_mol)
                # deep_print(mol)
                if a1.implicit_hydrogen_cnt() == 1:
                    h_atom = Atom("H", None)
                    atoms1 = atoms1[:1] + [h_atom] + atoms1[1:]
                    atoms2 = atoms2[:1] + [h_atom] + atoms2[1:]
                order_match = compare_stereochemistry(atoms1,atoms2)
                chirality_match = 1 * (sub_mol.chirality[a1] == mol.chirality[matches[a1]])
                if order_match != chirality_match:
                    return []
            #check allenal chirality:
            for c1,c2 in sub_mol.c_double_bonds:
                if c1 not in sub_mol.chirality or c2 not in sub_mol.chirality or matches[c1] not in mol.chirality or matches[c2] not in mol.chirality:
                    continue
                atoms_submol_c1 = [matches[n] for n,_ in sub_mol.ordered_bonds[c1]]
                atoms_mol_c1 = [n for n,_ in mol.ordered_bonds[matches[c1]]]
                atoms_submol_c2 = [matches[n] for n,_ in sub_mol.ordered_bonds[c2]]
                atoms_mol_c2 = [n for n,_ in mol.ordered_bonds[matches[c2]]]
                #treat implicit hydrogens as the second bonded atom:
                #http://opensmiles.org/spec/open-smiles-3-input.html
                if c1.implicit_hydrogen_cnt() == 1:
                    h_atom = Atom("H", None)
                    atoms_submol_c1 = atoms_submol_c1[:1] + [h_atom] + atoms_submol_c1[1:]
                    atoms_mol_c1 = atoms_mol_c1[:1] + [h_atom] + atoms_mol_c1[1:]
                if c2.implicit_hydrogen_cnt() == 1:
                    h_atom = Atom("H", None)
                    atoms_submol_c2 = atoms_submol_c2[:1] + [h_atom] + atoms_submol_c2[1:]
                    atoms_mol_c2 = atoms_mol_c2[:1] + [h_atom] + atoms_mol_c2[1:]
                order_match_c1 = compare_stereochemistry(atoms_submol_c1,atoms_mol_c1)
                chirality_match_c1 = 1 * (sub_mol.chirality[c1] == mol.chirality[matches[c1]])
                order_match_c2 = compare_stereochemistry(atoms_submol_c2,atoms_mol_c2)
                chirality_match_c2 = 1 * (sub_mol.chirality[c2] == mol.chirality[matches[c2]])
                if (order_match_c1 == chirality_match_c1) != (order_match_c2 == chirality_match_c2):
                    return []
            for neighbor, is_same_component_boolean in sub_mol.component_bonds.get(a1,[]):
                observed_same_component = (mol.atom2component[matches[a1]] == mol.atom2component[matches[neighbor]])
                if observed_same_component != is_same_component_boolean:
                    return []
        return [matches.copy()]

    #recursion
    out = []
    for a1 in sorted(sub_mol.explicit_atoms, key = lambda x: len(possibilities.get(x,[]))):
        if a1 in matches:
            continue
        if a1 not in possibilities:
            return []
        for m in possibilities[a1]:
            if m is not None and m in matches.values(): #don't map two atoms from sub_mol to the same atom in mol
                continue
            matches[a1] = m
            
            poss_cnt = 1
            for v in possibilities.values():
                poss_cnt *= len(v)
            if poss_cnt > 10000: #reprocess possibilities only if we have a lot because the function call is expensive
                next_poss = gen_possibilities(sub_mol,mol,matches)
            else:
                next_poss = possibilities
            out += match_substructure_helper(sub_mol, mol, next_poss, matches, done)
            matches.pop(a1)
        break #only process one atom from sub_mol.explicit_atoms, then break
    return out
########## end substructure matching ##########

class Search(object):
    @classmethod
    def search(self, end, start=None, verbose=False):
        q = Queue.PriorityQueue()
        i = (self.score(end, start),end,[])
        if verbose:
            print "Starting: ",(i[0], str(i[1]), i[2])
        done_list = set([]) #list of mols we have already checked
        q.put(i)
        while not q.empty():
            val,m,path = q.get()
            #check if we're already processed this reactant
            if m in done_list:
                continue
            done_list.add(m)
            if verbose:
                print done_list
                print "Searching... ", (val, str(m), path)
            if self.score(m,start) < -0.99:
                return val,str(m),list(reversed(path))
                break #synthesis complete
            for rxn in all_retros.values():
                if verbose:
                    print "Trying {rxn}...".format(**vars())
                out = rxn.run([m])
                if out:
                    for product_set in out:
                        for new_mol in product_set:
                            if new_mol:
                                if verbose:
                                    print new_mol
                                    print "Reaction result: " + str([str(new_mol)])
                                    print self.score(new_mol,start)
                                i = (self.score(new_mol,start) + 0.5 * len(path),new_mol,path+[str(rxn)])
                                q.put(i)
    @classmethod
    def score(self, end, start=None):
        """complexity score of a compound
        Should be high when that compound is difficult to synthesize
        and low when it is easy to synthesize
        """
        sim = end.similarity(start)
        if start:
            return -1 * sim #lower scores are better

        if Hydrobromination.has_br(end):
            return 2
        else:
            return 1

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
    print s, s_out
    #want to make sure smiles prints
    # CC(Br)CC instead the equivalent but
    # uglier CC(CC)Br
    assert(s_out == s)

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
    #mcmurray 6.8-6.9
    #TODO: use markovnikov constraint
    all_retros["hydrobromination"] = Retro("Hydrobromination", [Reaction('[C:1][C:2]([C:3])(Br)[C:4]>>[C:1][C:2]([C:3])=[C:4]'),Reaction('[C:1][C:2](Br)[C:3]([!C:4])([!C:5])>>[C:1][C:2]=[C:3]([*:4])([*:5])')])

    #mcmurray 7.1
    all_retros["dehydrohalogenation"] = Retro("Dehydrohalogenation", [Reaction('[*:1][C@:2]([*:3])=[C@@:4]([*:5])([*:6]) >> [*:1][C@@H:2]([*:3])[C@@:4]([*:5])([*:6])Br')])

    #mcmurray 7.1
    all_retros["alcohol_dehydration"] = Retro("Alcohol Dehydration", [Reaction('[C:1]=[C:2] >> [C:1][C:2]O')])

    #mcmurray 7.2
    all_retros["alkene_plus_x2"] = Retro("Alkene Plus X2", [Reaction('[*:2][C@@:1]([*:3])(Br)[C@:4](Br)([*:6])[*:5] >> [*:2][C@@:1]([*:3])=[C@@:4]([*:5])[*:6]')])
    
    #mcmurray 7.3
    # all_retros["halohydrin_formation"] = Retro("Halohydrin Formation", ['[C:1][C:2]([C:3])(Br)[C:4][OH]>>[C:1][C:2]([C:3])=[C:4]', '[C:1][C:2](Br)[C:3]([!C:4])([!C:5])[OH]>>[C:1][C:2]=[C:3]([*:4])([*:5])'])

    #mcmurray 7.3
    all_retros["halohydrin_formation"] = Retro("Halohydrin Formation", [Reaction('[*:1][C@@:2]([*:3])(Br)[C@:4]([*:6])([*:7])[OH]>>[*:1][C@@:2]([*:3])=[C@:4]([*:6])([*:7])', constraints = [Markovnikov("2","4")])])

    #mcmurray 7.4
    all_retros["oxymercuration"] = Retro("Oxymercuration", [Reaction('[C:1]([OH])[C:2]>>[C:1]=[C:2]', constraints = [Markovnikov("1","2")])])
    #mcmurray 7.5
    all_retros["hydroboration"] = Retro("Hydroboration", [Reaction('[*:1][C@@:3]([*:2])([H])[C@@:4]([*:5])([*:6])[O] >> [*:1][C@@:3]([*:2])=[C@:4]([*:5])[*:6]', constraints = [Markovnikov("3","4")])])

    #mcmurray 7.6
    all_retros["alkene_hydrogenation"] = Retro("Alkene Hydrogenation", [Reaction('[*:1][C@@:3]([*:2])([H])[C@@:4]([*:5])([*:6])[H] >> [*:1][C@@:3]([*:2])=[C@:4]([*:5])[*:6]')])

    #mucmurray 7.8
    all_retros["alkene_hydroxylation"] = Retro("Alkene Hydroxylation", [Reaction('[*:1][C@@:3]([*:2])([O])[C@@:4]([*:5])([*:6])[O] >> [*:1][C@@:3]([*:2])=[C@:4]([*:5])[*:6]')])

    #mcmurray 7.6
    all_retros["dichlorocarbene_addition"] = Retro("Dichlorocarbene Addition", [Reaction('[C:1]1C(Cl)(Cl)[C:2]1>>[C:1]=[C:2]')])

    #mcmurray 7.6
    all_retros["simmons_smith"] = Retro("Simmons-Smith", [Reaction('[C:1]1[CH2][C:2]1>>[C:1]=[C:2]')])

    #mcmurray 7.8
    all_retros["ozonolysis"] = Retro("Ozonolysis", [Reaction("[C:1]=O.[C:2]=O >> [C:1]=[C:2]")])

    # all_retros["esterification"] = Retro("esterification",["C(=O)O.OCC>>C(=O)OCC.O"]))

    #mcmurray 7.8
    all_retros["kmno4"] = Retro("KMnO4", [Reaction('[C:1]([*:2])([*:3])=O.[C:4]([*:5])([*:6])=O>>[C:1]([*:2])([*:3])=[C:4]([*:5])([*:6])'),
                                          Reaction('[C:3][C:1](=O)[OH].[C:2](=O)(=O) >> [C:3][CH:1]=[CH2:2]')])

    #mcmurray 7.9
    all_retros["diol_1_2_oxidative_cleavage"] = Retro("Diol 1 2 Oxidative Cleavage", [Reaction('[C:1]=O.[C:2]=O>>[C:1]([OH])[C:2]([OH])')])

    #mcmurray 8.3
    all_retros["dehydrohalogenation_vicinal_dihalides"] = Retro("Dehydrohalogenation Vicinal Dihalides", [Reaction('[C:1]#[C:2] >> Br[C:1][C:2]Br')])

    #mcmurray 8.9
    all_retros["acetylide_ion_alkylation"] = Retro("Acetylide Ion Alkylation", [Reaction('[CH:1]#[C:2][C:3][C:4] >> [CH:1]#[CH:2].[C:4][CH2:3][Br]'),
                                                                                Reaction('[C:1][C:2]#[C:3][CH2:4][C:5] >> [C:1][C:2]#[CH:3].[C:5][CH2:4][Br]')])

    #mcmurray 8.4
    all_retros["alkyne_plus_HX"] = Retro("Alkyne Plus HX", [Reaction('[C:1][C:2](Br)(Br)[CH3:3] >> [C:1][C:2](Br)=[CH2:3]'), Reaction('[C:1][C:2](Br)=[CH2:3] >> [C:1][C:2]#[CH:3]')])

    #mcmurray 8.4
    all_retros["alkyne_plus_X2"] = Retro("Alkyne Plus X2", [Reaction('[C:1][C:2](Br)(Br)[C:3](Br)(Br)[C:4] >> Br[C@@:2]([C:1])=[C@@:3]([C:4])Br'), Reaction('Br[C@@:2]([C:1])=[C@@:3]([C:4])Br >> [C:1][C:2]#[C:3][C:4]')])

    #mcmurray 8.5
    all_retros["alkyne_mercuric_hydration"] = Retro("Alkyne Mercuric Hydration", [Reaction('[C:1][C:2](=O)[CH3:3] >> [C:1][C:2]#[CH:3]')])

    #mcmurray 8.5
    all_retros["alkyne_hydroboration"] = Retro("Alkyne Hydroboration", [Reaction('[C:1][CH2:2][CH:3](=O) >> [C:1][C:2]#[CH:3]')])

    #mcmurray 8.6
    all_retros["alkyne_hydrogenation_paladium"] = Retro("Alkyne Hydrogenation Paladium", [Reaction('[C:1][CH2:2][CH2:3][C:4] >> [C:1][C:2]#[C:3][C:4]')])

    #mcmurray 8.6
    all_retros["alkyne_hydrogenation_lindlar"] = Retro("Alkyne Hydrogenation Lindlar", [Reaction('[C:1][C@@H:2]=[C@@H:3][C:4] >> [C:1][C:2]#[C:3][C:4]')])

    #mcmurray 8.6
    all_retros["alkyne_hydrogenation_lithium"] = Retro("Alkyne Hydrogenation Lithium", [Reaction('[C:1][C@@H:2]=[C@H:3][C:4] >> [C:1][C:2]#[C:3][C:4]')])

    #mcmurray 8.8
    #create acetylide anions?
    
    #mcmurray 8.8
    all_retros["acetylide_ion_alkylation"] = Retro("Acetylide Ion Alkylation", [Reaction('[CH:1]#[C:2][CH2:3][C:4] >> [CH:1]#[CH:2].[C:4][CH2:3]Br'), Reaction('[C:1][C:2]#[C:3][CH2:4][C:5]  >> [C:1][C:2]#[CH:3].[C:5][CH2:4]Br')])

    #mcmurray 8.7
    # all_retros["alkyne_oxidative_cleavage"] = Retro("Alkyne Oxidative Cleavage", [Reaction('[*:1][C:2](=O)[OH].[*:4][C:3](=O)[OH] >> [*:1][C:2]#[C:3][*:4]')])
    all_retros["alkyne_oxidative_cleavage"] = Retro("Alkyne Oxidative Cleavage", [Reaction('[*:2][C:1](=O)[OH] >> C#[C:1][*:2]')])

    #test reaction that just inverts the chirality of a chiral carbon
    all_retros["invert_chirality"] = Retro("Invert Chirality", [Reaction('[*:1][C@]([*:2])([*:3])[*:4] >> [*:1][C@@]([*:2])([*:3])[*:4]')])
    
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
    assert(all_retros["hydrobromination"].run([start1])[0] == (Molecule("CC(C)=C"),))

    #dehydrohalogenation
    start1 = (Molecule("I[C@@](Cl)=[C@@](I)Cl"),)
    end1 = (Molecule("I[C@@](Cl)[C@@](Cl)(I)Br"),)
    assert(all_retros["dehydrohalogenation"].run(start1)[0] == end1)

    
    #alcohol_dehydration
    start1 = (Molecule("C=CCCC"),)
    end = ((Molecule("CCCCCO"),),(Molecule("CCCC(O)C"),))
    out = all_retros["alcohol_dehydration"].run(start1)
    assert(set(end) == set(out))

    #[*:2][C@@:1](F)([*:3])[C@:4](F)([*:6])[*:5] >> [*:2][C:1]([*:3])=[C:4]([*:5])[*:6]
    start1 = (Molecule("C[C@@](Cl)(Br)[C@](C)(Cl)Br"),)
    end1 = (Molecule("Cl[C@](C)=[C@](C)Cl"),)
    assert(all_retros["alkene_plus_x2"].run(start1)[0] == end1)
    
    start1 = (Molecule("CC(C)(Br)CO"),)
    end1 = (Molecule("CC(C)=C"),)
    assert(all_retros["halohydrin_formation"].run(start1)[0] == end1)

    start1 = (Molecule("CC(O)C(Br)(I)"),)
    end1 = (Molecule("C(C)=C(Br)I"),)
    assert(all_retros["oxymercuration"].run(start1)[0] == end1)

    start1 = (Molecule("C(Br)(Br)C"),)
    end1 = (Molecule("C(Br)(Br)=C"),)
    assert(all_retros["alkene_hydrogenation"].run(start1)[0] == end1)

    start1 = (Molecule("BrC(Br)(Br)C(Br)(Br)"),)
    end1 = start1
    assert(all_retros["alkene_hydroxylation"].run(start1)[0] == end1)
    
    start1 = (Molecule("OCCO"),)
    end1 = (Molecule("C=C"),)
    assert(all_retros["alkene_hydroxylation"].run(start1)[0] == end1)

    start1 = (Molecule("CC(Cl)(Cl)C"),)
    end1 = (Molecule("CC(Cl)(Cl)C"),)
    assert(all_retros["dichlorocarbene_addition"].run(start1)[0] == end1)

    start1 = (Molecule("C1C(Cl)(Cl)C1"),)
    end1 = (Molecule("C=C"),)
    assert(all_retros["dichlorocarbene_addition"].run(start1)[0] == end1)
    
    start1 = (Molecule("CCC1CC1"),)
    end1 = (Molecule("CCC=C"),)
    assert(all_retros["simmons_smith"].run(start1)[0] == end1)

    start1 = [Molecule("BrC=O"),Molecule("ClC=O")]
    end1 = (Molecule("ClC=CBr"),)
    assert(all_retros["ozonolysis"].run(start1)[0] == end1)

    start = (Molecule("CC=CC(C)C"),)
    end = (Molecule("CC=O.CC(C)C=O"),)
    assert(all_retros["ozonolysis"].run(end)[0][0] == start[0])
    
    start1 = [Molecule("C(C)(C)=O"),Molecule("C(C)(C)=O")]
    end1 = (Molecule("CC(C)=C(C)C"),)
    assert(all_retros["kmno4"].run(start1)[0] == end1)

    start1 = [Molecule("C=O"),Molecule("C=O")]
    end1 = (Molecule("C(O)C(O)"),)
    assert(all_retros["diol_1_2_oxidative_cleavage"].run(start1)[0] == end1)

    start1 = (Molecule("CC#CC"),)
    end1 = (Molecule("CC(Br)C(Br)C"),)
    assert(all_retros["dehydrohalogenation_vicinal_dihalides"].run(start1)[0] == end1)

    start1 = (Molecule("CC#CCC"),)
    end1 = (Molecule("CC#C"),Molecule("CCBr"))
    out1 = all_retros["acetylide_ion_alkylation"].run(start1)[0]
    assert(set(out1) == set(end1))

    start1 = (Molecule("CC(Br)(Br)C"),)
    end1 = (Molecule("CC(Br)=C"),)
    assert(all_retros["alkyne_plus_HX"].run(start1)[0] == end1)

    start1 = (Molecule("Br[C@@](C)=[C@@](C)Br"),)
    end1 = (Molecule("CC#CC"),)
    assert(all_retros["alkyne_plus_X2"].run(start1)[0] == end1)

    start1 = (Molecule("CC(=O)C"),)
    end1 = (Molecule("CC#C"),)
    assert(all_retros["alkyne_mercuric_hydration"].run(start1)[0] == end1)

    start1 = (Molecule("CCC=O"),)
    end1 = (Molecule("CC#C"),)
    assert(all_retros["alkyne_hydroboration"].run(start1)[0] == end1)

    start1 = (Molecule("CCCC"),)
    end1 = (Molecule("CC#CC"),)
    assert(all_retros["alkyne_hydrogenation_paladium"].run(start1)[0] == end1)
    
    start1 = (Molecule("C[C@@]=[C@@]C"),)
    end1 = (Molecule("CC#CC"),)
    assert(all_retros["alkyne_hydrogenation_lindlar"].run(start1)[0] == end1)

    start1 = (Molecule("C[C@@]=[C@]C"),)
    end1 = (Molecule("CC#CC"),)
    assert(all_retros["alkyne_hydrogenation_lithium"].run(start1)[0] == end1)

    start1 = (Molecule("C#CCC"),)
    end1 = (Molecule("C#C"),Molecule("CCBr"))
    out1 = all_retros["acetylide_ion_alkylation"].run(start1)[0]
    assert(set(out1) == set(end1))

    #alkyne_oxidative_cleavage simplified for now -- R group is just a hydrogen
    start1 = (Molecule("CC(=O)O"),)
    end1 = (Molecule("CC#C"),)
    assert(all_retros["alkyne_oxidative_cleavage"].run(start1)[0] == end1)

    
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

def test_inter_intra_component():
    s = Substructure("(CCCC.CC)")
    m = Molecule("CCCC.CC")
    assert(not find_substructure(s,m))

    s = Substructure("CCCC.CC")
    m = Molecule("CCCC.CC")
    assert(find_substructure(s,m))

    s = Substructure("(CCCC).(CC)")
    m = Molecule("CCCC.CC")
    assert(find_substructure(s,m))

    s = Substructure("(CCCC).(CC)")
    m = Molecule("CCCCCC")
    assert(not find_substructure(s,m))

    s = Substructure("CCCC.CC")
    m = Molecule("CCCCCC")
    assert(find_substructure(s,m))


def test_merge_split():
    m1 = Molecule("CCC")
    m2 = Molecule("CCC")
    m3 = Molecule("CCC")
    mout = Molecule("CCC.CCC.CCC")
    assert(Molecule.merge([m1,m2,m3]) == mout) #yikes this is slow

    m1 = Molecule("C1CC1")
    m2 = Molecule("C2CC2")
    m3 = Molecule("C3CC3")
    mout = Molecule("C1CC1.C2CC2.C3CC3")
    assert(Molecule.merge([m1,m2,m3]) == mout) #yikes this is slow

    mout = Molecule("C2CC2.CCCI.CC(Br)(Br)C")
    assert(set(mout.split()) == set([Molecule("C2CC2"),Molecule("CCCI"),Molecule("CC(Br)(Br)C")]))

def test_cis_trans():
    m1 = Molecule('I[C@@](Cl)=[C@](F)Br')
    m2 = Molecule('I[C@](Cl)=[C@@](F)Br')
    m3 = Molecule('I[C@@](Cl)=[C@@](F)Br')
    m4 = Molecule('I[C@](Cl)=[C@](F)Br')
    assert(m1 == m2)
    assert(m1 != m3)
    assert(m1 != m4)
    assert(m2 != m3)
    assert(m2 != m4)
    assert(m3 == m4)
    
def test_search():
    #McMurray 7.26 (pg 239)
    
    #Alkene Hydroxylation (X)
    start = Molecule('C1CC=CC1')
    end = Molecule('C1C([OH])C([OH])CC1')
    score, m_out, path = Search.search(end, start)
    print(score, m_out, path)
    assert(score == -1 and path == ['Alkene Hydroxylation'])
    
    #Hydroboration (currently blocked by markovnikov condition)
    # start = Molecule('C1CC=CC1')
    # end = Molecule('C1CC([OH])CC1')

    #Dichlorocarbene Addition (X)
    start = Molecule('C1CC=CC1')
    end = Molecule('C1CC2C(Cl)(Cl)C2C1')

    score, m_out, path = Search.search(end, start)
    assert(score == -1 and path == ["Dichlorocarbene Addition"])
    
    #Alcohol dehydration (X)
    start = Molecule('C1C([CH3])([OH])CCCC1')
    end = Molecule('C1C([CH3])=CCCC1')
    score, m_out, path = Search.search(end, start)
    assert(score == -1 and path == ["Alcohol Dehydration"])

    #Hydroboration (X)
    start = Molecule('CC(C)=C')
    end = Molecule('CC(C)C[OH]')

    score, m_out, path = Search.search(end, start)
    assert(score == -1 and path == ["Hydroboration"])
    
    ####alkyne test cases
    #McMurray 8.29 (pg 271)

    #alkyne mercuric hydration (X)
    start = Molecule('[CH3][CH2]C#[CH]')
    end = Molecule('[CH3][CH2]C(=O)[CH3]')

    score, m_out, path = Search.search(end, start)
    assert(score == -1 and path == ["Alkyne Mercuric Hydration"])

    #alkyne hydroboration (X)
    start = Molecule('[CH3][CH2]C#[CH]')
    end = Molecule('CCCC=O')
    score, m_out, path = Search.search(end, start)
    assert(score == -1 and path == ["Alkyne Hydroboration"])

    start = Molecule('CCC#C')
    end = Molecule('CCC(O)=O')
    score, m_out, path = Search.search(end, start)
    assert(score == -1 and path == ["Alkyne Oxidative Cleavage"])

    start = Molecule('CCCCC=C')
    end = Molecule('CCCCC#C')
    score, m_out, path = Search.search(end, start)
    assert(path == ['Alkene Plus X2', 'Dehydrohalogenation Vicinal Dihalides'])

    
def test_implicit_hydrogen_stereochemistry():
    #implicit H is treated as the second neighbor of the carbon
    m = Molecule("Br[C@](Cl)I")
    m2 = Molecule("Br[C@@](Cl)I")
    assert(m != m2)

    m2 = Molecule("Br[C@@](I)Cl")
    assert(m == m2)

def deep_print(obj):
    def deep_print_helper(obj):
        """Recursively call str() on all sub-parts of an object
        example: a list of objects
        """
        if isinstance(obj,tuple):
            return tuple([deep_print_helper(i) for i in obj])
        elif isinstance(obj,dict):
            return "dict(" + str([(deep_print_helper(k),deep_print_helper(v)) for k,v in obj.items()]) + ")"
        elif isinstance(obj,list):
            return [deep_print_helper(i) for i in obj]
        elif isinstance(obj,set):
            return "set(" + str([deep_print_helper(i) for i in obj]) + ")"
        else:
            return str(obj)
    print deep_print_helper(obj)

    
if __name__ == "__main__":
    rxn_setup()

    # test_search()

    # start = Molecule("CC(C)CC#C")
    # end = Molecule("CC(C)")
    # score, m_out, path = Search.search(end, start, verbose=True)


    # start = Molecule('CCCCC=C')
    # end = Molecule('CCCCC#C')
    # score, m_out, path = Search.search(end, start, verbose=True)
    
    # start = Molecule("CCCC#C")
    # end = Molecule("CCCC=O")
    # score, m_out, path = Search.search(end, start, verbose=True)


    # start = Molecule("CCCCC#C")
    # end = Molecule("CCCC[C@]1C[C@@]1C")
    # score, m_out, path = Search.search(end, start, verbose=True)

    start = Molecule('CCCCC=C')
    end = Molecule('CCCCC#C')
    score, m_out, path = Search.search(end, start, verbose=True)
    assert(path == ['Alkene Plus X2', 'Dehydrohalogenation Vicinal Dihalides'])
    
