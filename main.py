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
import random

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


## user can enter keyword
keyword = "Glycolysis"
# keyword = "Citric acid cycle (TCA cycle)"
# keyword = "Signaling by Insulin receptor"
# keyword = "Cell Cycle"

def parse_keyword(keyword):

    # print("keyword is ", keyword)
    for i in range(len(keyword.split())):
        if i == 0:
            str_key = keyword.split()[i]
        else:
            str_key = str_key + "%20" + keyword.split()[i]
    return str_key

# this block of code helps us to get the ID for the keyword entered
def get_id_pathcommons(keyword):
    #print(f'https://reactome.org/ReactomeRESTfulAPI/RESTfulWS/biopaxExporter/Level3/{Id}')
    # res = requests.get('http://www.pathwaycommons.org/pc2/search.json?q='+keyword+'&datasource=reactome')
    keyword_changed = keyword.title()
    #print(parse_keyword(keyword_changed))
    url = 'http://www.pathwaycommons.org/pc2/search.json?q='+parse_keyword(keyword_changed)+'&datasource=reactome'
    print(url)
    res = requests.get(url)
    # print(res.text)
    if res.ok:
        return res.text
        # this returns a long python string that is in JSON format
    else:
        print(res.ok)
        return None


def extract_id_from_search(id_info_text, keyword):
    for i in range(len(json.loads(id_info_text)["searchHit"])):
        if keyword in str(json.loads(id_info_text)["searchHit"][i]['name']):
            print(str(json.loads(id_info_text)["searchHit"][i]['name']))
            id_url = json.loads(id_info_text)["searchHit"][i]['uri']
            break
        else:
            print("cannot find")
    num = ""
    for c in id_url:
        if c.isdigit():
            num = num + c
    id = int(num)
    return id
names = dict()



### UX functions
@eel.expose
def extract_id_from_search2(keyword):

    global names
    names = dict()
    id_info_text = get_id_pathcommons(keyword)
    for i in range(len(json.loads(id_info_text)["searchHit"])):
        if keyword.lower() in str(json.loads(id_info_text)["searchHit"][i]['name']).lower():
            print(str(json.loads(id_info_text)["searchHit"][i]['name']))
            id_url = json.loads(id_info_text)["searchHit"][i]['uri']
            names[str(json.loads(id_info_text)["searchHit"][i]['name'])] = id_url

    return list(names.keys())

@eel.expose
def generate_text_id(id):
    print(f"working on query with id {id}")
    pathway_dict = {}
    search = get_id_pathcommons(id)
    id = extract_id_from_search(search, id)
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
        reaction_buzzword = row['buzzwords']

        pathway_dict[row['biochemicalReactionDispName']] = create_descriptions(reaction_name, reaction_type, cellular_loc, reactant_list, product_list, enzyme_name, reaction_buzzword)
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


def find_catalysis_control(reaction):
  # search through Catalysis first for enzymes/facilitators. If nothing is found, also search through Control for facilitators.
  all_cat = list()
  for cat in catalysis:
    if cat.controlled == reaction:
      all_cat.append(cat)
  if len(all_cat) == 0:
    for con in control:
      if con.controlled == reaction:
        all_cat.append(con)
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

def determine_event_type(reaction):
  ## parses reaction.display_name AND reaction.comment string for keywords to determine event type.
  ## assigns reaction 3 scores. whichever category of score is the highest, reaction gets assigned that type.
  ## also returns a list of "buzzwords" to be used in NLG paragraph

  # 1. metabolic = catabolism AND anabolism (creation/production AND destruction of reactions & products)
  metabolic_symbols = ['=>', "<=>", "<=","+"]
  metabolic_keywords = ['anabolic', 'anabolism', 'biosynthesis', 'catabolic', 'catabolism', 
                          'catalyse', 'catalysed', 'catalyses', 'catalysing', 'catalysis', 
                          'catalyze', 'catalyzed', 'catalyzes', 'catalyzing', 'cleavage', 
                          'cleave', 'cleaved', 'cleaves', 'cleaving', 'conversion', 
                          'convert', 'converted', 'converting', 'converts', 'dehydrogenation', 'dehydrogenates',
                          'enzymatically', 'enzyme', 'enzymes', 'form', 'formation', 
                          'formed', 'forming', 'forms', 'generate', 'generated', 'generates', 
                          'generating', 'generation', 'hydrolyse', 'hydrolysed', 'hydrolyses', 'hydrolysing',
                          'hydrolysis', 'hydrolyze', 'hydrolyzed', 'hydrolyzes', 'hydrolyzing', 
                          'interconversion', 'intermediate', 'intermediates', 'intraconversion', 
                          'isomerisation', 'isomerises', 'metabolic', 'metabolism', 'oxidate', 
                          'oxidated', 'oxidates', 'oxidation', 'oxidating', 'oxidatively', 'oxidize', 'oxidizes', 'oxidized', 'oxidizing', 
                          'produce', 'produced', 'produces', 'producing', 'production', 'react', 'reacted', 
                          'reacting', 'reaction', 'reactions', 'reacts', 'reduce', 'reduces', 'reduced',
                          'reducing', 'reduction', 'reductively', 'synthesis', 'synthesize', 
                          'synthesized', 'synthesizes', 'synthesizing', 'yield', 'yielded', 'yielding', 'yields']

  metabolic_buzzwords = {
      'anabolism': ['anabolic', 'anabolism'],
      'biosynthesis': ['biosynthesis'],
      'catabolism': ['catabolic', 'catabolism'],
      'catalysis': ['catalyse', 'catalysed', 'catalyses', 'catalysing', 'catalysis', 'catalyze', 'catalyzed', 'catalyzes', 'catalyzing'],
      'cleavage': ['cleavage', 'cleave', 'cleaved', 'cleaves', 'cleaving'],
      'conversion': ['conversion', 'convert', 'converted', 'converting', 'converts'],
      'dehydrogenation': ['dehydrogenation', 'dehydrogenates'],
      'formation': ['form', 'formation', 'formed', 'forming', 'forms'],
      'generation': ['generate', 'generated', 'generates', 'generating', 'generation'],
      'hydrolysis': ['hydrolyse', 'hydrolysed', 'hydrolyses', 'hydrolysing', 'hydrolysis', 'hydrolyze', 'hydrolyzed', 'hydrolyzes', 'hydrolyzing'],
      'interconversion': ['interconversion'],
      'intraconversion': ['intraconversion'],
      'isomerisation': ['isomerisation', 'isomerises'],
      'metabolism': ['metabolic', 'metabolism'],
      'oxidation': ['oxidate', 'oxidated', 'oxidates', 'oxidation', 'oxidating', 'oxidatively', 'oxidize', 'oxidizes', 'oxidized', 'oxidizing'],
      'production': ['produce', 'produced', 'produces', 'producing', 'production'],
      'reduction': ['reduce', 'reduces', 'reduced', 'reducing', 'reduction', 'reductively'],
      'synthesis': ['synthesis', 'synthesize', 'synthesized', 'synthesizes', 'synthesizing']
  }

  # 2. signaling = activation AND deactivation. also includes regulation/regulatory events. a shift from one state to another.
  signaling_keywords = ['activate', 'activated', 'activates', 'activating', 'activation', 'activator', 
                          'activators', 'active', 'activities', 'activity', 'agonist', 'agonists', 
                          'antagonist', 'antagonists', 'assemble', 'assembled', 'assembles', 'assembling', 
                          'assembly', 'associate', 'associated', 'associates', 'associating', 'association', 
                          'bind', 'binding', 'binds', 'bound', 'cascade', 'cascades', 'cascading', 'conform', 
                          'conformation', 'conforming', 'conforms', 'control', 'controlled', 'controller', 
                          'controlling', 'controls', 'de-phosphorylate', 'de-phosphorylated', 'de-phosphorylates', 
                          'de-phosphorylating', 'de-phosphorylation', 'deactivate', 'deactivated', 'deactivates', 
                          'deactivating', 'deactivation', 'degradation', 'degrade', 'degraded', 'degrades', 'degrading', 
                          'dephosphorylate', 'dephosphorylated', 'dephosphorylates', 'dephosphorylating', 'dephosphorylation', 
                          'dissociate', 'dissociated', 'dissociates', 'dissociating', 'dissociation', 'downstream', 
                          'express', 'expressed', 'expresses', 'expressing', 'expression', 'glycosylate', 'glycosylated', 
                          'glycosylates', 'glycosylating', 'glycosylation', 'inactivate', 'inactivated', 'inactivates', 
                          'inactivating', 'inactivation', 'inactive', 'inhibit', 'inhibited', 'inhibiting', 'inhibition', 
                          'inhibitor', 'inhibitors', 'inhibitory', 'inhibits', 'initiate', 'initiated', 'initiates', 'initiating', 
                          'initiation', 'initiator', 'initiators', 'interact', 'interacting', 'interaction', 'interactor', 
                          'interactors', 'interacts', 'ligand', 'ligands', 'mediate', 'mediated', 'mediates', 'mediating', 
                          'mediation', 'modification', 'modifications', 'modified', 'modifies', 'modify', 'modifying', 'phosphorylate', 
                          'phosphorylated', 'phosphorylates', 'phosphorylating', 'phosphorylation', 'receive', 'received', 'receives', 
                          'receiving', 'receptor', 'receptors', 'regulate', 'regulated', 'regulates', 'regulating', 'regulation', 
                          'regulator', 'regulatory', 'repress', 'repressed', 'represses', 'repressing', 'repression', 'signal', 'signaled', 
                          'signaling', 'signalled', 'signalling', 'signals', 'stimulate', 'stimulated', 'stimulates', 'stimulating', 
                          'stimulation', 'transcribe', 'transcribed', 'transcribes', 'transcribing', 'transcription', 'translate', 
                          'translated', 'translates', 'translating', 'translation', 'trigger', 'triggered', 'triggering', 'triggers', 
                          'ubiquitin', 'ubiquitinate', 'ubiquitinated', 'ubiquitinates', 'ubiquitinating', 'ubiquitination', 'upstream',
                          'deubiquitination', 'deubiquitinate', 'deubiquitinates', 'deubiquitinated', 'deubiquitinating']

  signaling_buzzwords = {
      'activation': ['activate', 'activated', 'activates', 'activating', 'activation', 'activator', 'activators', 
                      'active', 'activities', 'activity', 'agonist', 'agonists',
                      'initiate', 'initiated', 'initiates', 'initiating', 'initiation', 'initiator', 'initiators',
                      'stimulate', 'stimulated', 'stimulates', 'stimulating', 'stimulation'],
      'assembly': ['assemble', 'assembled', 'assembles', 'assembling', 'assembly'],
      'association': ['associate', 'associated', 'associates', 'associating', 'association'],
      'binding': ['bind', 'binding', 'binds', 'bound', 'ligand', 'ligands',
                  'receive', 'received', 'receives', 'receiving', 'receptor', 'receptors'],
      'signalling cascade': ['cascade', 'cascades', 'cascading','trigger', 'triggered', 'triggering', 'triggers'],
      'degradation': ['degradation', 'degrade', 'degraded', 'degrades', 'degrading'],
      'dephosphorylation': ['de-phosphorylate', 'de-phosphorylated', 'de-phosphorylates', 'de-phosphorylating', 'de-phosphorylation',
                              'dephosphorylate', 'dephosphorylated', 'dephosphorylates', 'dephosphorylating', 'dephosphorylation'],
      'dissociation': ['dissociate', 'dissociated', 'dissociates', 'dissociating', 'dissociation'],
      'inhibition': ['antagonist', 'antagonists', 'deactivate', 'deactivated', 'deactivates', 'deactivating', 'deactivation',
                          'inactivate', 'inactivated', 'inactivates', 'inactivating', 'inactivation', 'inactive', 'inhibit', 'inhibited', 'inhibiting', 'inhibition', 
                          'inhibitor', 'inhibitors', 'inhibitory', 'inhibits'],
      'gene expression': ['express', 'expressed', 'expresses', 'expressing', 'expression', 'transcribe', 'transcribed', 'transcribes', 'transcribing', 'transcription', 'translate', 
                          'translated', 'translates', 'translating', 'translation'],
      'gene repression': ['repress', 'repressed', 'represses', 'repressing', 'repression'],
      'glycosylation': ['glycosylate', 'glycosylated', 'glycosylates', 'glycosylating', 'glycosylation'],
      'modification': ['modification', 'modifications', 'modified', 'modifies', 'modify', 'modifying'],
      'phosphorylation': ['phosphorylate', 'phosphorylated', 'phosphorylates', 'phosphorylating', 'phosphorylation'],
      'regulation': ['control', 'controlled', 'controller', 'controlling', 'controls', 'downstream', 'mediate', 'mediated', 'mediates', 'mediating', 'mediation',
                      'regulate', 'regulated', 'regulates', 'regulating', 'regulation', 'regulator', 'regulatory', 'upstream'],
      'ubiquitination': ['ubiquitin', 'ubiquitinate', 'ubiquitinated', 'ubiquitinates', 'ubiquitinating', 'ubiquitination'],
      'deubiquitination': ['deubiquitination', 'deubiquitinate', 'deubiquitinates', 'deubiquitinated', 'deubiquitinating']
  }

  # 3. transport = movement of something from one place to another
  transport_keywords = ['endocytose', 'endocytosed', 'endocytoses', 'endocytosing', 'endocytosis', 
                          'enter', 'entered', 'entering', 'enters', 'exit', 'exited', 'exiting', 'exits', 
                          'exocytose', 'exocytosed', 'exocytoses', 'exocytosing', 'exocytosis', 'export', 
                          'exported', 'exporting', 'exports', 'externalisation', 'externalise', 'externalised', 
                          'externalises', 'externalising', 'externalization', 'externalize', 'externalized', 
                          'externalizes', 'externalizing', 'import', 'imported', 'importing', 'imports', 'internalisation', 
                          'internalise', 'internalised', 'internalises', 'internalising', 'internalization', 'internalize', 
                          'internalized', 'internalizes', 'internalizing', 'leave', 'leaves', 'leaving', 'localization', 
                          'localize', 'localized', 'localizes', 'localizing', 'secrete', 'secreted', 'secretes', 
                          'secreting', 'secretion', 'secretory', 'traffic', 'trafficked', 'trafficking', 'traffics', 
                          'translocate', 'translocated', 'translocates', 'translocating', 'translocation', 
                          'transport', 'transportation', 'transported', 'transporting', 'transports']

  transport_buzzwords = {
      'endocytosis': ['endocytose', 'endocytosed', 'endocytoses', 'endocytosing', 'endocytosis'],
      'exocytosis': ['exocytose', 'exocytosed', 'exocytoses', 'exocytosing', 'exocytosis'],
      'export': ['export', 'exported', 'exporting', 'exports'],
      'externalization': ['externalisation', 'externalise', 'externalised', 'externalises', 'externalising', 'externalization', 'externalize', 'externalized', 'externalizes', 'externalizing'],
      'import': ['import', 'imported', 'importing', 'imports'],
      'internalization': ['internalisation', 'internalise', 'internalised', 'internalises', 'internalising', 'internalization', 'internalize', 'internalized', 'internalizes', 'internalizing'],
      'localization': ['localization', 'localize', 'localized', 'localizes', 'localizing'],
      'secretion': ['secrete', 'secreted', 'secretes', 'secreting', 'secretion', 'secretory'],
      'translocation': ['translocate', 'translocated', 'translocates', 'translocating', 'translocation']
  }
  
  # 4. miscellaneous event --> assign as signaling to be safe

  def get_buzz(word, buzzwords_dict):
    buzz = "NONE"
    for k in buzzwords_dict:
        for i in range(len(buzzwords_dict[k])):
            if buzzwords_dict[k][i] == word:
                buzz = k
    return buzz

  # actual score tallying
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
      if w in metabolic_symbols:
          score_metabolic += 1
      if w in metabolic_keywords:
          score_metabolic += 2
      if w in signaling_keywords:
          score_signaling += 2
      if w in transport_keywords:
          score_transport += 3

  met = ("metabolic", score_metabolic, metabolic_buzzwords)
  sig = ("signaling", score_signaling, signaling_buzzwords)
  tran = ("transport", score_transport, transport_buzzwords)

  if met[1] == sig[1] == tran[1] == 0:
    type = sig[0]
    buzz_dict = sig[2]
  else:
    winner = max(met[1],sig[1],tran[1])
    if winner is met[1]:
      type = met[0]
      buzz_dict = met[2]
    elif winner is tran[1]:
      type = tran[0]
      buzz_dict = tran[2]
    else:
      type = sig[0]
      buzz_dict = sig[2]

  # accumulates list of buzzwords to be returned
  buzzwords = []

  for r in reaction.display_name.lower().split():
    buzz = get_buzz(r, buzz_dict)
    if buzz != "NONE" and buzz not in buzzwords:
      buzzwords.append(buzz)
  for w in reaction.comment[0].lower().split():
    buzz = get_buzz(w, buzz_dict)
    if buzz != "NONE" and buzz not in buzzwords:
      buzzwords.append(buzz)

  return (type, buzzwords)


### Generate Data Frame based on pybiopax model
def generate_df(model):

    reactions = get_reactions(model)
    catalysis = get_instances_list(model, 'catalysis')
    control = get_instances_list(model, 'control')
    df = pd.DataFrame(
        columns=["biochemicalReactionDispName", "reactionType", "cellularLoc", "reactants", "products", "enzyme","buzzwords"])

    if len(reactions) == 0:
        print("no reactions")
        return df

    for reaction in reactions:
        cats = find_catalysis_control(reaction)
        listEnzymes = list()
        listEnzymeNames = list()
        listEnzymeLoc = list()
        if len(cats) > 0:
            for cat in cats:
                listEnzymes += get_protein(cat)
        if len(listEnzymes) > 0:
            for enz in listEnzymes:
                listEnzymeNames.append(get_prot_name(enz))
                if get_loc(enz).term[0] not in listEnzymeLoc:
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
        (reaction_Type, reaction_Buzz) = determine_event_type(reaction)

        dic = {
            "biochemicalReactionDispName": reaction.display_name,
            "reactionType": reaction_Type,
            "cellularLoc": listEnzymeLoc,
            "reactants": listLeft,
            "products": listRight,
            "enzyme": listEnzymeNames,
            "buzzwords": reaction_Buzz
        }
        df.loc[len(df.index)] = dic
    return df


#df = generate_df(model)


##### Natural Language Generation #####

def gen_sentence_list(entity_list):
    if len(entity_list) == 1:
        gen_str = entity_list[0]
    elif len(entity_list) > 2:
        gen_str = entity_list[0]
        for ent in entity_list[1:-1]:
            gen_str = gen_str + ", " + ent
        gen_str = gen_str + ", and " + entity_list[-1]
    else:
        gen_str = entity_list[0]+ " and " + entity_list[-1]
    return gen_str


def create_descriptions(reaction_name, reaction_type, cell_loc, reactant_list, product_list, enzyme, buzzwords):
    para = ''

    # sentence 1
    sentence1 = nlgFactory.createClause()
    sentence1_noun = nlgFactory.createNounPhrase("This")

    if reaction_type == "metabolic":
        reaction_type_noun = nlgFactory.createNounPhrase("a reaction")
        reac_procduct_verb = nlgFactory.createVerbPhrase("has")
        reac_proc_part1 = nlgFactory.createNounPhrase("the reactants")
        reac_proc_part1.addPostModifier(gen_sentence_list(reactant_list))
        reac_proc_part2 = nlgFactory.createNounPhrase("the products")
        reac_proc_part2.addPostModifier(gen_sentence_list(product_list))

    elif reaction_type == "transport":
        reaction_type_noun = nlgFactory.createNounPhrase("a transport event")

        sub_for_involve = ["involves", "consists of", "includes"]
        reac_procduct_verb = nlgFactory.createVerbPhrase(random.choice(sub_for_involve))

        reac_proc_part1 = nlgFactory.createNounPhrase("the transport of")
        if reactant_list == product_list:
            reac_proc_part1.addPostModifier(gen_sentence_list(reactant_list))
        else:
            reac_proc_part1.addPostModifier(
                gen_sentence_list(reactant_list) + ", and also " + gen_sentence_list(product_list))
    else:
        reaction_type_noun = nlgFactory.createNounPhrase("an event")
        # reac_procduct_verb = nlgFactory.createVerbPhrase("involves")

        sub_for_involve = ["involves", "consists of", "includes"]
        reac_procduct_verb = nlgFactory.createVerbPhrase(random.choice(sub_for_involve))

        reac_proc_part1 = nlgFactory.createNounPhrase(gen_sentence_list(reactant_list))
        reac_proc_part1.addPreModifier("the entities")
        reac_proc_part1.addPostModifier("in its initial state")
        reac_proc_part2 = nlgFactory.createNounPhrase(gen_sentence_list(product_list))
        reac_proc_part2.addPostModifier("the entities")
        reac_proc_part2.addComplement("in its secondary state")

    reaction_type_noun.addPostModifier("that")
    sub_place_in = ["happens in", "takes place in", "occurs in"]
    cellloc_verb = nlgFactory.createVerbPhrase(random.choice(sub_place_in))
    # cellloc_verb.addPostModifier('place in')

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