#!/usr/bin/env python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import DataStructs
import rdkit.Chem.rdchem

import Queue
import sys
import time

"""
Unsolved rdkit problems:
(1) Stereochemistry doesn't work for reactions: the reactant template matches regardless of stereochemistry and the stereochemistry information is only used for the difference in chirality between the reactant and the product.

(2) Hard to specify Atom:1 == Atom:2
For example want a pattern to match ClCCCl, BrCCBr, ICCI, but not ClCCBr, etc.

(3) Can't specify how substituted a carbon is
In particular you can match carbon with C, can match non-carbon, non-hydrogen with !C but there's no way to match all non-carbon atoms (ie including hydrogen!)

(4) (minor) Multiple reactants have to be in the right order for the reaction to run to completion
"""

def plot_mol(m):
    """needed to install imagemagick on ubuntu for this to work"""
    img = Draw.MolToImage(m)
    img.show()

def show_smiles(s):
    m = Chem.MolFromSmiles(s)
    plot_mol(m)

def score(end, start=None):
    """complexity score of a compound
    Should be high when that compound is difficult to synthesize
    and low when it is easy to synthesize
    """
    f_end = FingerprintMols.FingerprintMol(end)
    f_start = FingerprintMols.FingerprintMol(start)
    sim = DataStructs.FingerprintSimilarity(f_end, f_start)
    if start:
        return -1 * sim #lower scores are better
    
    if Hydrobromination.has_br(end):
        return 2
    else:
        return 1

class RXN(object):
    def __init__(self):
        pass
    def reverse(self, m):
        pass


# class Electrophilic_Addition_HX(RXN):

def run_rxn_list(rxn_list, reactants):
    """reactants is a tuple"""
    if not isinstance(reactants,tuple):
        raise
    products = []
    for rxn in rxn_list:
        products += rxn.RunReactants(reactants)
    if products:
        return products
    else:
        return (reactants,)

class Reaction(object):
    def __init__(self, name, rxn_string_list):
        self.__name__ = name
        self.rxn_list = [AllChem.ReactionFromSmarts(r) for r in rxn_string_list]
    def __str__(self):
        return self.name
    def reverse(self, reactants):
        return Reaction.run_once(reactants,self.rxn_list)
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
    def run_once(cls, reactants, rxn_list):
        try:
            products = [p for product_set in run_rxn_list(rxn_list, tuple(reactants)) for p in product_set]
        except ValueError, e:
            #invalid number of reactants
            return reactants
        for m in products:
            print Chem.MolToSmiles(m)
            Chem.SanitizeMol(m)
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
            new_print = FingerprintMols.FingerprintMol(m)
            if any([new_print == p for p in prints]):
                continue
            prints.append(new_print)
            out.append(m)
        return out
    @classmethod
    def equal(cls, mol_list1, mol_list2):
        return len(cls.drop_duplicates(mol_list1+ mol_list2)) == len(mol_list1)

    
class Hydrobromination(Reaction):
    @classmethod
    def has_br(self, m):
        patt = Chem.MolFromSmarts('Br')
        match = m.GetSubstructMatch(patt)
        return match
    @classmethod
    def forward(self, m):
        pass
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1][C:2]([C:3])(Br)[C:4]>>[C:1][C:2]([C:3])=[C:4]')
        rxn2 = AllChem.ReactionFromSmarts('[C:1][C:2](Br)[C:3]([!C:4])([!C:5])>>[C:1][C:2]=[C:3]([*:4])([*:5])')
        products = Reaction.run_once(m,[rxn1,rxn2])
        return products


class Dehydrohalogenation(Reaction):
    """CC(Br) --> C=C"""
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1]=[C:2] >> [C:1][C:2][F,Cl,Br,I]')
        products = Reaction.run_once(m, [rxn1])
        return products

class AlcoholDehydration(Reaction):
    """CCOH --> C=C"""
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1]=[C:2] >> [C:1][C:2][O;H1]')
        products = Reaction.run_once(m, [rxn1])
        return products

class Alkene_Plus_X2(Reaction):
    """C=C -> XCCX"""
    @classmethod
    def reverse(self, m):
        rxn_pattern = '[*:2][C@@:1](F)([*:3])[C@:4](F)([*:6])[*:5]'
        if len(m) == 1 and m[0].HasSubstructMatch(Chem.MolFromSmarts(rxn_pattern), useChirality=True, recursionPossible=True): #no match in ring
            rxn1 = AllChem.ReactionFromSmarts('{rxn_pattern} >> [*:2][C:1]([*:3])=[C:4]([*:5])[*:6]'.format(**vars()))
            products = Reaction.run_once(m,[rxn1])
            return products
        return m
        
    @classmethod
    def reverse_old(self,m):    
        #DANGER: can't figure out how to make ReactionFromSmarts match stereochemistry properly.
        #It's currently matching regardless of the stereochemistry of the input.
        #So I'm just filtering using HasSubstructMatch for now, though that'll fail when for
        #example, both reactant stereoisomers exist in the same molecule

        ring_match = '[C;R](F)(*)[C;R](F)(*)'
        ring_rxn_pattern = '[*;R:2][C@@:1](F)([*:3])[C@:4](F)([*:6])[*;R:5]'
        nonring_rxn_pattern = '[*:2][C@@:1](F)([*:3])[C@:4](F)([*:6])[*:5]'
        patt = Chem.MolFromSmarts(reactant_smarts)
        #recursionPossible=True allows matching the 'substructure' as the whole molecule
        #http://www.rdkit.org/Python_Docs/rdkit.Chem.rdchem.Mol-class.html#HasSubstructMatch
        if m.HasSubstructMatch(Chem.MolFromSmarts(ring_match), useChirality=True, recursionPossible=True):
            if m.HasSubstructMatch(Chem.MolFromSmarts(ring_rxn_pattern), useChirality=True, recursionPossible=True):
                rxn1 = AllChem.ReactionFromSmarts('{ring_rxn_pattern} >> [*:2]/[C:1]([*:3])=[C:4]([*:5])\[*:6]'.format(**vars()))
                products = Reaction.run_once([m],[rxn1])
                return products
        elif m.HasSubstructMatch(Chem.MolFromSmarts(nonring_rxn_pattern), useChirality=True, recursionPossible=True): #no match in ring
            rxn1 = AllChem.ReactionFromSmarts('{nonring_rxn_pattern} >> [*:2]/[C:1]([*:3])=[C:4]([*:5])\[*:6]'.format(**vars()))
            products = Reaction.run_once([m],[rxn1])
            return products

        return m

class Halohydrin_Formation(Reaction):
    """C=C -> XCCOH"""
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1][C:2]([C:3])(Br)[C:4][OH]>>[C:1][C:2]([C:3])=[C:4]')
        rxn2 = AllChem.ReactionFromSmarts('[C:1][C:2](Br)[C:3]([!C:4])([!C:5])[OH]>>[C:1][C:2]=[C:3]([*:4])([*:5])')
        return Reaction.run_once(m,[rxn1, rxn2])

class Oxymercuration(Reaction):
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1][C:2]([C:3])([OH])[CH:4]>>[C:1][C:2]([C:3])=[C:4]')
        rxn2 = AllChem.ReactionFromSmarts('[C:1][C:2]([OH])[CH:3]([!C:4])([!C:5])>>[C:1][C:2]=[C:3]([*:4])([*:5])')
        return Reaction.run_once(m,[rxn1, rxn2])


class Hydroboration(Reaction):
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1][C:2]([C:3])[CH:4]([OH])>>[C:1][C:2]([C:3])=[C:4]')
        rxn2 = AllChem.ReactionFromSmarts('[C:1][C:2][CH2:3][OH]>>[C:1][C:2]=[C:3]')
        rxn3 = AllChem.ReactionFromSmarts('[C:1][CH:2]([OH])[CH2:3][C:4]>>[C:1][C:2]=[C:3][C:4]')
        return Reaction.run_once(m,[rxn1,rxn2, rxn3])

class Alkene_Hydrogenation(Reaction):
    @classmethod
    def reverse(self, m):
        #for some reason the commented reaction crashes on 'CCC(C)=O'
        #so instead using rxn1-5 for now
        # rxn = AllChem.ReactionFromSmarts('[C;!H0:1][C;!H0:2]>>[C:1]=[C:2]')
        rxn1 = AllChem.ReactionFromSmarts('[CH2:1][CH3:2]>>[CH1:1]=[CH2:2]')
        rxn2 = AllChem.ReactionFromSmarts('[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]')
        rxn3 = AllChem.ReactionFromSmarts('[CH2:1][CH2:2]>>[CH:1]=[CH:2]')
        rxn4 = AllChem.ReactionFromSmarts('[CH2:1][CH:2]>>[CH:1]=[CH0:2]')
        rxn5 = AllChem.ReactionFromSmarts('[CH:1][CH:2]>>[CH0:1]=[CH0:2]')
        return Reaction.run_once(m,[rxn1,rxn2,rxn3,rxn4,rxn5])

class Alkene_Hydroxylation(Reaction):
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1]([OH])[C:2]([OH])>>[C:1]=[C:2]')
        return Reaction.run_once(m,[rxn1])

class DichlorocarbeneAddition(Reaction):
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1]C(Cl)(Cl)[C:2]>>[C:1]=[C:2]')
        return Reaction.run_once(m,[rxn1])

class SimmonsSmith(Reaction):
    @classmethod
    def reverse(self, m):
        rxn1 = AllChem.ReactionFromSmarts('[C:1]1[CH2][C:2]1>>[C:1]=[C:2]')
        return Reaction.run_once(m,[rxn1])

class Ozonolysis(Reaction):
    @classmethod
    def reverse(self, reactants):
        rxn1 = AllChem.ReactionFromSmarts('[C:1]=O.[C:2]=O >> [C:1]=[C:2]')
        return Reaction.run_once(reactants,[rxn1])

class KMnO4(Reaction):
    @classmethod
    def reverse(self, reactants):
        rxn1 = AllChem.ReactionFromSmarts('[C:1]([C:2])([C:3])=O.[C:4]([C:5])([C:6])=O>>[C:1]([C:2])([C:3])=[C:4]([C:5])([C:6])')
        rxn2 = AllChem.ReactionFromSmarts('[C:3][C:1](=O)[OH].[C:2](=O)(=O) >> [C:3][CH:1]=[CH2:2]')
        return Reaction.run_once(reactants,[rxn1, rxn2])

class Diol_1_2_Oxidative_Cleavage(Reaction):
    @classmethod
    def reverse(self, reactants):
        rxn1 = AllChem.ReactionFromSmarts('[C:1]=O.[C:2]=O >> [C:1]([OH])[C:2]([OH])')
        return Reaction.run_once(reactants, [rxn1])


    
RXN_LIST = [Hydrobromination, Dehydrohalogenation, AlcoholDehydration, Alkene_Plus_X2, Halohydrin_Formation, Oxymercuration, Hydroboration, Alkene_Hydrogenation, Alkene_Hydroxylation, DichlorocarbeneAddition, SimmonsSmith, Ozonolysis, KMnO4, Diol_1_2_Oxidative_Cleavage,
            Reaction('Dehydrohalogenation_Vicinal_Dihalides', ['[C:1][C:2]#[C:3][C:4] >> [C:1][C:2](Br)[C:3](Br)[C:4]',
                                                               '[C:1][C:2]#[C:3][C:4] >> [C:1][C:2]=[C:3](Br)[C:4]']),
            Reaction('Acetylide_Ion_Alkylation',['[CH:1]#[C:2][C:3][C:4] >> [CH:1]#[CH:2].[C:4][CH2:3][Br]',
                                                 '[C:1][C:2]#[C:3][CH2:4][C:5] >> [C:1][C:2]#[CH:3].[C:5][CH2:4][Br]']),
            Reaction('Alkyne_Plus_HX', ['[C:1][C:2](Br)(Br)[CH3:3] >> [C:1][C:2](Br)=[CH2:3]',
                                        '[C:1][C:2](Br)=[CH2:3] >> [C:1][C:2]#[CH:3]']),
            Reaction('Alkyne_Plus_X2', ['[C:1][C:2](Br)(Br)[C:3](Br)(Br)[C:4] >> Br\[C:2]([C:1])=[C:3]([C:4])/Br',
                                        'Br\[C:2]([C:1])=[C:3]([C:4])/Br >> [C:1][C:2]#[C:3][C:4]']),
            Reaction('Alkyne_Mercuric_Hydration', ['[C:1][C:2](=O)[CH3:3] >> [C:1][C:2]#[CH:3]']),
            Reaction('Alkyne_Hydroboration', ['[C:1][CH2:2][CH:3](=O) >> [C:1][C:2]#[CH:3]']),
            Reaction('Alkyne_Hydrogenation_Paladium', ['[C:1][CH2:2][CH2:3][C:4] >> [C:1][C:2]#[C:3][C:4]']),
            Reaction('Alkyne_Hydrogenation_Lindlar', ['[C:1]\[CH:2]=[CH:3]/[C:4] >> [C:1][C:2]#[C:3][C:4]']),
            Reaction('Alkyne_Hydrogenation_Lithium', ['[C:1]\[CH:2]=[CH:3]\[C:4] >> [C:1][C:2]#[C:3][C:4]']),
            Reaction('Alkyne_To_Acetylene', ['[C:1][C:2]#[C-:3] >> [C:1][C:2]#[CH:3]']),
            Reaction('Acetylide_Ion_Alkylation',['[CH:1]#[C:2][CH2:3][C:4] >> [CH:1]#[CH:2].[C:4][CH2:3]Br',
                                                 '[C:1][C:2]#[C:3][CH2:4][C:5]  >> [C:1][C:2]#[CH:3].[C:5][CH2:4]Br']),
            Reaction('Alkyne_Oxidative_Cleavage', ['[C:1][C:2](=O)[OH].[C:4][C:3](=O)[OH] >> [C:1][C:2]#[C:3][C:4]'])
]


            
# RXN_LIST = [Hydroboration]

class Search(object):
    @classmethod
    def search(self, end, start=None):
        q = Queue.PriorityQueue()
        i = (score(end, start),end,[])
        print "Starting: ",(i[0], Chem.MolToSmiles(i[1]), i[2])
        done_list = [] #list of mols we have already checked
        q.put(i)
        while not q.empty():
            val,m,path = q.get()
            skip = False
            #check if we're already processed this reactant
            for d in done_list:
                if FingerprintMols.FingerprintMol(m) == FingerprintMols.FingerprintMol(d):
                    skip = True #already processed this row
                    break
            if skip:
                continue
            done_list.append(m)
            print done_list
            print "Searching... ", (val, Chem.MolToSmiles(m), path)
            # show_mol(m)
            # raw_input("Waiting for input...")
            if score(m,start) < -0.99:
                return val,Chem.MolToSmiles(m),list(reversed(path))
                break #synthesis complete
            for rxn in RXN_LIST:
                print "Trying {rxn.__name__}...".format(**vars())
                out = rxn.reverse([m])
                if out:
                    for new_mol in out:
                        if new_mol:
                            print "Reaction result: " + str([Chem.MolToSmiles(new_mol)])
                            print score(new_mol,start)
                            i = (score(new_mol,start),new_mol,path+[rxn.__name__])
                            q.put(i)

def basic_test():
    m = Chem.MolFromSmiles('c1nccc2n1ccc2')
    AllChem.Compute2DCoords(m)
    template = Chem.MolFromSmiles('c1nccc2n1ccc2')
    AllChem.Compute2DCoords(template)
    AllChem.GenerateDepictionMatching2DStructure(m,template)
    m = Chem.MolFromSmiles("O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5")


def test_search():
    m = Chem.MolFromSmiles("CCCCC(C)C(Br)CCCCCCBr")
    score, m_out, path =Search.search(m)
    print "Synthesis solved!"
    print "Starting molecule:"
    print Chem.MolToSmiles(m_out)
    print "Sequence of reactions:"
    print path
    # show_mol(m_out)


def test(rxn_string, reactant_string_list):
    reactants = [Chem.MolFromSmiles(s) for s in reactant_string_list]
    for r in reactants:
        plot_mol(r)
    time.sleep(2)
    fn = globals()[rxn_string]
    print [Chem.MolToSmiles(p) for p in fn.reverse(reactants)]
    print [plot_mol(p) for p in fn.reverse(reactants)]

def readCL():
    parser = argparse.ArgumentParser()
    parser.add_argument("-l","--log_level",default="INFO")
    args = parser.parse_args()
    return args.log_level
    
if __name__ == "__main__":
    m = Chem.MolFromSmiles('CC[C@@H]([OH])[C@]([Br])(I)I')
    # m = Chem.MolFromSmiles("CCCCC(C)C(Br)C([OH])CCCCCBr")
    m = Chem.MolFromSmiles('[OH]CC(Br)(Br)')

    m = Chem.MolFromSmiles('CC1C(Cl)(Cl)C1C')

    m1 = Chem.MolFromSmiles('CCCC(=O)[OH]')
    m2 = Chem.MolFromSmiles('O=C=O')


    m1 = Chem.MolFromSmiles('CC(=O)C')
    m2 = Chem.MolFromSmiles('CC(=O)C')

    
    start = Chem.MolFromSmiles('C1CC=CC1')
    end = Chem.MolFromSmiles('C1C([OH])C([OH])CC1')

    start = Chem.MolFromSmiles('C1CC=CC1')
    end = Chem.MolFromSmiles('C1CC([OH])CC1')

    start = Chem.MolFromSmiles('C1CC=CC1')
    end = Chem.MolFromSmiles('C1CC2C(Cl)(Cl)C2C1')

    start = Chem.MolFromSmiles('C1C([CH3])([OH])CCCC1')
    end = Chem.MolFromSmiles('C1C([CH3])=CCCC1')

    # start = Chem.MolFromSmiles('CC=CC(C)C')
    # end = 

    start = Chem.MolFromSmiles('CC(C)=C')
    end = Chem.MolFromSmiles('CC(C)C[OH]')


    ####alkyne test cases
    start = Chem.MolFromSmiles('[CH3][CH2]C#[CH]')
    end = Chem.MolFromSmiles('[CH3][CH2]C(=O)[CH3]')


    start = Chem.MolFromSmiles('[CH3][CH2]C#[CH]')
    end = Chem.MolFromSmiles('CCCC=O')


    start = Chem.MolFromSmiles('CCCC#C')
    end = Chem.MolFromSmiles('CCCC=O')
    
    print Search.search(end, start)
