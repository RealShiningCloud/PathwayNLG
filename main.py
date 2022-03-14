### KR 550 Final Project: NLG with Reactions (using Reactome and BioPax)
## Serena Jinchen Xie, Ehsan Alipour, Sharon Wong
## March 2022

##### Import Dependencies #####
# pip3 install requests, pybiopax, pandas, simplenlg, eel

import requests
import json
import pybiopax
import pandas as pd
import simplenlg

from pybiopax.paths import find_objects
from pybiopax.biopax.model import get_sub_objects
from pybiopax.biopax.interaction import Catalysis
from pybiopax.biopax.interaction import Control

from simplenlg.framework import *
from simplenlg.lexicon import *
from simplenlg.realiser.english import *
from simplenlg.phrasespec import *
from simplenlg.features import *

import eel

lexicon = Lexicon.getDefaultLexicon()
nlgFactory = NLGFactory(lexicon)
realiser = Realiser(lexicon)



##### Editable Variables #####

# list of reactions / id:
Id_gly = 70171  # glycolysis (metabolic)
Id_cca = 71403  # citric acid cycle (metabolic)

Id_wnt = 195721  # WNT signaling pathways (signaling)
Id_ins = 74752  # insulin (signaling)
Id_egfr = 177929  # EGFR signaling pathways (signaling)

Id_trans = 392499  # metabolism of proteins? (includes dna transcription / translation) (challenge)
Id_dna = 69306  # dna replication (challenge)
Id_cell = 1640170  # cell cycle (challenge)
model = ""

##### Query Reactome API -> return Reaction ID# from Keyword #####

id = 1640170


### UX functions
@eel.expose
def generate_text_id(id):
    print(f"working on query with id {id}")
    pathway_dict = {}
    reactome = get_reactome(id)
    if reactome == None:
        print("ID not found.")
        return "Failed 1"

    global model
    print("Biopax Model Loaded")
    df = generate_df(model)
    if len(df.index) == 0:
        print("Error loading reactions")
        return "Failed 2"
    df.reset_index()
    print("Generating Text...")
    for index, row in df.iterrows():
        reaction_name = row['biochemicalReactionDispName']
        reaction_type = row['reactionType']
        cellular_loc = row['cellularLoc']
        reactant_list =row['reactants']
        product_list =row['products']
        enzyme_name = row['enzyme']

        pathway_dict[row['biochemicalReactionDispName']] = create_descriptions(reaction_name, reaction_type, cellular_loc, reactant_list, product_list, enzyme_name)
        listNames = list(pathway_dict.keys())
        listTexts = list(pathway_dict.values())
    print("Done.")
    return listTexts

@eel.expose
def get_reaction_names():
    global model
    listNames = list()

    reactions = get_reactions(model)
    for re in reactions:
        listNames.append(re.display_name)
    return listNames

@eel.expose
def get_title():
    return[model.objects['Pathway1'].display_name, model.objects['Pathway1'].comment]
### Extracting contents thru Reactome API

def get_reactome(Id):
    # print(f'https://reactome.org/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/{Id}')
    print("Sending Query...")
    res = requests.get(f'https://reactome.org/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/{Id}')
    global model
    # print(res.text)
    if res.ok:
        print("Query Successful!")
        file = open('biopax.owl', 'w')
        file.write(res.text)
        file.close()
        model = pybiopax.model_from_owl_file("biopax.owl")
        get_all_data(model)
        assert isinstance(model.objects, dict)
        assert all(isinstance(obj, pybiopax.biopax.BioPaxObject)
                   for obj in model.objects.values())
        return res.text
    else:
        print("Query failed.")
        print(res.ok)
        return None


##### Create pybiopax Model & Extract Elements #####

# documentation: https://github.com/indralab/pybiopax
'''
model = pybiopax.model_from_owl_file("biopax.owl")
assert isinstance(model.objects, dict)
assert all(isinstance(obj, pybiopax.biopax.BioPaxObject)
           for obj in model.objects.values())

pathway_name = model.objects['Pathway1'].display_name
pathway_descrip = model.objects['Pathway1'].comment[0]
'''

##### pybiopax Helper Functions #####

def get_reactions(model):
    # generally, return only the Reactions in Pathway1. if no Reactions in Pathway1, return ALL reactions.
    reactions = list()
    for entity in model.objects['Pathway1'].pathway_component:
        if isinstance(entity, pybiopax.biopax.interaction.BiochemicalReaction):
            reactions.append(entity)
    if len(reactions) == 0:
        for entity_key in model.objects:
            if isinstance(model.objects[entity_key], pybiopax.biopax.interaction.BiochemicalReaction):
                reactions.append(model.objects[entity_key])
    return reactions


def get_instances_list(model, instanceStr):
    all_reactions = list()
    for instance in model.objects.keys():
        if instanceStr in instance.lower():
            all_reactions.append(model.objects[instance])
    return all_reactions


def find_catalysis(reaction):
    all_cat = list()
    for cat in catalysis:
        if cat.controlled == reaction:
            all_cat.append(cat)
    return all_cat


def get_protein(cat, model=model):
    return list(cat.controller)


def get_loc(prot, model=model):
    return prot.cellular_location


def get_prot_name(prot, model=model):
    return prot.display_name


##### END pybiopax Helper Functions #####
def get_all_data(model):
    global reactions
    global catalysis
    global control
    reactions = get_reactions(model)
    catalysis = get_instances_list(model, 'catalysis')
    control = get_instances_list(model, 'control')
    return True


#### Parse and Determine Type #####

def determine_event_type(reaction):
    ## parses reaction.display_name AND reaction.comment string for keywords to determine event type.
    ## assigns reaction 3 scores. whichever category of score is the highest, reaction gets assigned that type.

    # 1. metabolic = catabolism AND anabolism (creation/production AND destruction of reactions & products)
    metabolic_symbols = ['=>', "<=>", "<=", "+"]
    metabolic_keywords = ["catalysis", "catalyzes", "catalyze", "catalyzing",
                          "catalyses", "catalyse", "catalysing",
                          "synthesis", "synthesizes", "synthesize", "synthesizing", "biosynthesis", "intermediate",
                          "intermediates",
                          "conversion", "converts", "convert", "converting", "interconversion", "intraconversion",
                          "formation", "forms", "form", "forming", "generation", "generates", "generate", "generating",
                          "production", "produces", "produce", "producing",
                          "hydrolysis", "hydrolysed", "hydrolyses", "isomerises", "isomerisation",
                          "oxidate", "oxidates", "oxidation", "oxidizing", "oxidatively", "reduce", "reduces",
                          "reduction", "reducing", "reductively",
                          "cleavage", "cleaves", "cleave", "cleaving",
                          "catabolic", "catabolism", "anabolic", "anabolism"]

    # 2. signaling = activation AND deactivation. also includes regulation/regulatory events. a shift from one state to another.
    signaling_keywords = ["signal", "signaling", "signals", "signaled", "ligand", "ligands",
                          "receptor", "receptors", "receive", "receives", "received", "receiving",
                          "regulatory", "regulator", "regulates", "regulate", "regulating", "regulated",
                          "bind", "binds", "binding", "bound", "conform", "conforms", "conformation", "conforming",
                          "dissociate", "dissociates", "dissociation", "associate", "associates", "association",
                          "inactive", "active", "inhibit", "inhibits", "inhibiting", "inhibitor", "inhibitory",
                          "activate", "activates", "activating", "activator",
                          "inactivation", "activation", "repression", "represses", "repressed", "repressing",
                          "control", "controls", "controlled", "controlling", "controller",
                          "agonist", "antagonist", "interact", "interacts", "interaction", "interactor", "interactors",
                          "phosphorylation", "phosphorylate", "phosphorylates", "phosphorylating", "phosphorylated",
                          "dephosphorylation", "dephosphorylate", "dephosphorylates", "dephosphorylating",
                          "dephosphorylated",
                          "glycosylation", "glycosylated", "glycosylate", "glycolsylates",
                          "expression", "expressed", "expressing", "transcription", "translation", "translated",
                          "transcribed",
                          "assembly", "assembled", "assembles", "degradation", "degraded", "degrading", "degrades",
                          "ubiquitination", "ubiquitinated", "ubiquitinates", "ubiquitinate", "ubiquitinating"]

    # 3. transport = movement of something from one place to another
    transport_keywords = ["transport", "transports", "transportation", "transported", "transporting",
                          "translocate", "translocates", "translocation", "translocated", "translocating",
                          "export", "import", "exported", "imported", "exports", "imports", "leaves", "enters",
                          "localization", "localize", "localizes", "localizes",
                          "exocytosis", "endocytosis", "exocytosed", "endocytosed",
                          "traffic", "traffics", "trafficked", "trafficking",
                          "secretion", "secreted", "secretes", "secrete", "secretory",
                          "internalization", "internalized", "internalize", "internalizes"]

    # 4. miscellaneous event --> assign as signaling to be safe

    score_metabolic = 0
    score_signaling = 0
    score_transport = 0

    for r in reaction.display_name.lower().split():
        if r in metabolic_symbols:
            score_metabolic += 2
        if r in metabolic_keywords:
            score_metabolic += 3
        if r in signaling_keywords:
            score_signaling += 4
        if r in transport_keywords:
            score_transport += 5

    for w in reaction.comment[0].lower().split():
        if r in metabolic_symbols:
            score_metabolic += 1
        if r in metabolic_keywords:
            score_metabolic += 2
        if r in signaling_keywords:
            score_signaling += 2
        if r in transport_keywords:
            score_transport += 3

    met = ("metabolic", score_metabolic)
    sig = ("signaling", score_signaling)
    tran = ("transport", score_transport)

    if met[1] == sig[1] == tran[1] == 0:
        type = sig[0]
    else:
        winner = max(met[1], sig[1], tran[1])
        if winner is met[1]:
            type = met[0]
        elif winner is tran[1]:
            type = tran[0]
        else:
            type = sig[0]

    return type


### Generate Data Frame based on pybiopax model
def generate_df(model):

    reactions = get_reactions(model)
    catalysis = get_instances_list(model, 'catalysis')
    control = get_instances_list(model, 'control')
    df = pd.DataFrame(
        columns=["biochemicalReactionDispName", "reactionType", "cellularLoc", "reactants", "products", "enzyme"])

    if len(reactions) == 0:
        print("no reactions")
        return df

    for reaction in reactions:
        cats = find_catalysis(reaction)
        listEnzymes = list()
        listEnzymeNames = list()
        listEnzymeLoc = list()
        if len(cats) > 0:
            for cat in cats:
                listEnzymes += get_protein(cat)
        if len(listEnzymes) > 0:
            for enz in listEnzymes:
                listEnzymeNames.append(get_prot_name(enz))
                listEnzymeLoc.append(get_loc(enz).term[0])

        # get full names of left and right (reactants & products) as a list of string
        listLeft = list()
        listRight = list()
        for i in range(len(reaction.left)):
            if len(reaction.left[i].name) > 1:
                listLeft.append(reaction.left[i].name[1])
            else:
                listLeft.append(reaction.left[i].name[0])

        for j in range(len(reaction.right)):
            if len(reaction.right[j].name) > 1:
                listRight.append(reaction.right[j].name[1])
            else:
                listRight.append(reaction.right[j].name[0])

        # parse reaction name / display name to determine event type
        reaction_Type = determine_event_type(reaction)

        dic = {
            "biochemicalReactionDispName": reaction.display_name,
            "reactionType": reaction_Type,
            "cellularLoc": listEnzymeLoc,
            "reactants": listLeft,
            "products": listRight,
            "enzyme": listEnzymeNames,
        }
        df.loc[len(df.index)] = dic
    return df


#df = generate_df(model)


##### Natural Language Generation #####

# hi I can't get the newline character to work properly someone please help
# I also don't know if this function is necessary
def create_introduction(pathwayName, pathwayDescrip):
    title = "Pathway: " + pathwayName + '\n' + "Description: " + pathwayDescrip
    return title


def gen_sentence_list(reactant_list):
    if len(reactant_list) == 1:
        reactant_type = str(type(reactant_list[0])).split('.')[3][:-2]
        gen_str = reactant_type + ' ' + reactant_list[0].display_name
    elif len(reactant_list) > 2:
        gen_str = str(type(reactant_list[0])).split('.')[3][:-2] + ' ' + reactant_list[0].display_name
        for r in reactant_list[1:-1]:
            gen_str = gen_str + ", " + str(type(r)).split('.')[3][:-2] + ' ' + r.display_name
        gen_str = gen_str + ", and " + str(type(reactant_list[-1])).split('.')[3][:-2] + ' ' + reactant_list[
            -1].display_name
    else:
        gen_str = str(type(reactant_list[0])).split('.')[3][:-2] + ' ' + reactant_list[0].display_name
        gen_str = gen_str + " and " + str(type(reactant_list[1])).split('.')[3][:-2] + ' ' + reactant_list[
            -1].display_name

    return gen_str


def create_descriptions(reaction_name, reaction_type, cell_loc, reactant_list, product_list, enzyme):
    para = "Name: " + reaction_name

    # sentence 1
    sentence1 = nlgFactory.createClause()
    sentence1_noun = nlgFactory.createNounPhrase("This")

    if reaction_type == "metabolic":
        reaction_type_noun = nlgFactory.createNounPhrase("a reaction")
        reac_procduct_verb = nlgFactory.createVerbPhrase("has")
        reac_proc_part1 = nlgFactory.createNounPhrase("the reactants")
        reac_proc_part1.addPostModifier(gen_react_list(reactants_list))
        reac_proc_part2 = nlgFactory.createNounPhrase("the products")
        reac_proc_part2.addPostModifier(gen_react_list(products_list))
    elif reaction_type == "transport":
        reaction_type_noun = nlgFactory.createNounPhrase("a transport event")
        reac_procduct_verb = nlgFactory.createVerbPhrase("involves")
        reac_proc_part1 = nlgFactory.createNounPhrase("the transport of")
        if reactant_list == product_list:
            reac_proc_part1.addPostModifier(gen_sentence_list(reactant_list))
        else:
            reac_proc_part1.addPostModifier(
                gen_sentence_list(reactant_list) + ", and also " + gen_sentence_list(product_list))
    else:
        reaction_type_noun = nlgFactory.createNounPhrase("an event")
        reac_procduct_verb = nlgFactory.createVerbPhrase("involves")
        reac_proc_part1 = nlgFactory.createNounPhrase(gen_react_list(reactants_list))
        reac_proc_part1.addPreModifier("the entities")
        reac_proc_part1.addPostModifier("in its initial state")
        reac_proc_part2 = nlgFactory.createNounPhrase(gen_react_list(products_list))
        reac_proc_part2.addPostModifier("the entities")
        reac_proc_part2.addComplement("in its secondary state")

    reaction_type_noun.addPostModifier("that")
    cellloc_verb = nlgFactory.createVerbPhrase("takes")
    cellloc_verb.addPostModifier('place in')

    if len(cell_loc) != 0:
        cellloc_clause = nlgFactory.createClause()
        cellloc_noun = nlgFactory.createNounPhrase(cell_loc[0])
        cellloc_clause.setVerb(cellloc_verb)
        cellloc_clause.setObject(cellloc_noun)
        reaction_type_noun.addModifier(cellloc_clause)
        sentence1.setSubject(sentence1_noun)
        sentence1.setObject(reaction_type_noun)
        sentence1.setVerb("is")
        para = para + '\n' + realiser.realiseSentence(sentence1)

        # second sentence
        reac_procduct_sentence = nlgFactory.createClause()

        if reaction_type == "transport":
            obj_part = reac_proc_part1
        else:
            obj_part = nlgFactory.createCoordinatedPhrase(reac_proc_part1, reac_proc_part2)

        reac_procduct_sentence.setObject(obj_part)
        reac_procduct_sentence.setSubject('It')
        reac_procduct_sentence.setVerb(reac_procduct_verb)
        para = para + ' ' + realiser.realiseSentence(reac_procduct_sentence)

    elif len(cell_loc) == 0:
        reac_procduct_clause = nlgFactory.createClause()

        if reaction_type == "transport":
            obj_part = reac_proc_part1
        else:
            obj_part = nlgFactory.createCoordinatedPhrase(reac_proc_part1, reac_proc_part2)
        reac_procduct_clause.setObject(obj_part)

        reac_procduct_clause.setVerb(reac_procduct_verb)
        reaction_type_noun.addModifier(reac_procduct_clause)
        sentence1.setSubject(sentence1_noun)
        sentence1.setObject(reaction_type_noun)
        sentence1.setVerb("is")
        para = para + '\n' + realiser.realiseSentence(sentence1)

    if len(enzyme) == 0:
        print('here')
        last_sentence = "It's not catalyzed by any enzyme."
        para = para + ' ' + last_sentence
    else:
        last_sentence = nlgFactory.createClause()
        if reaction_type == "metabolic":
            last_sen_noun = nlgFactory.createNounPhrase("This reaction")
        elif reaction_type == 'transport':
            last_sen_noun = nlgFactory.createNounPhrase("This transport event")
        else:
            last_sen_noun = nlgFactory.createNounPhrase("This event")
        last_sen_verb = nlgFactory.createVerbPhrase("is facilitated by")
        last_sen_noun_2 = nlgFactory.createNounPhrase(enzyme[0])
        enzyme_name = nlgFactory.createNounPhrase("enzyme")

        last_sentence.setSubject(last_sen_noun)
        last_sentence.setObject(last_sen_noun)
        last_sentence.setVerb(last_sen_verb)
        last_sentence.setObject(last_sen_noun_2)
        if reaction_type == "metabolic":
            last_sentence.setIndirectObject(enzyme_name)

        para = para + ' ' + realiser.realiseSentence(last_sentence)
    return para


def main():
    eel.init('web')
    eel.start('index.html')

main()