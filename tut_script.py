import cobra.test
import warnings
import math

oxygen="o2"

acetate="ac"
acetaldehyde="acald"
two_oxoglutarate="akg"
ethanol="etoh"
d_fructose="fru"
fumarate="fum"
d_glucose="glc__D"
l_glutamine="gln__L"
l_glutamate="glu__L"
d_lactate="lac__D"
l_malate="mal__L"
pyruvate="pyr"
succinate="succ"
ATP="atp"
ADP="adp"
ATPM="ATPM"
proton="h"
threePG="3pg"

cytosolic_metabolite_suffix="_c"

glucose_exchange_realistic_lower_bound=-18.5
organic_exchange_realistic_lower_bound=-20.0
unlimited_reaction_bound=1000.0
limited_reaction_bound=0.0

class FBAExperiment:
    def __init__(self,name,reactions_bounds=[],aerobic=True,objective_function_name=None):
        self.name=name
        self.reactions_bounds=list(reactions_bounds)
        self.objective_function_name=objective_function_name
        if aerobic:
            self.name=self.name+"_aerobic"
            self.reactions_bounds.append(UnlimitedSubstrateReactionBounds(oxygen))
        else:
            self.name=self.name+"_unaerobic"
            self.reactions_bounds.append(UnavailableSubstrateReactionBounds(oxygen)) 
        
    def run(self,model):
        return self.set_model_bounds_and_get_optimal_values(model)
        
    def set_model_bounds_and_get_optimal_values(self,model):
        for reaction_bound in self.reactions_bounds:
            reaction_bound.set_reaction_bounds(model)
        if self.objective_function_name:
            model.objective=model.reactions.get_by_id(self.objective_function_name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            solution=model.optimize("max")
            return(solution)
      
    def run_experiment_and_print_result(self,model):
        pass
    
            
class SubstrateFBAExperiment(FBAExperiment):
    def __init__(self,substrate,aerobic):
        super(SubstrateFBAExperiment,self).__init__(name="Realistic_"+substrate,aerobic=aerobic)
        self.reactions_bounds.append(RealisticSubstrateReactionBounds(substrate))
        if not substrate==d_glucose:
            self.reactions_bounds.append(UnavailableSubstrateReactionBounds(d_glucose))
    
    def run(self,model):
        return super(SubstrateFBAExperiment,self).run(model)
        
    def run_experiment_and_print_result(self,model):
        opt_solution=self.run(model)
        print(self.name+": "+str(round(opt_solution.objective_value,4))+"/h")
        
class CofactorAndPrecursorsFBATest(FBAExperiment):
    drain_suffix="_drain"
    
    def __init__(self,target,aerobic):
        super(CofactorAndPrecursorsFBATest,self).__init__(name=target+CofactorAndPrecursorsFBATest.drain_suffix,aerobic=aerobic)
        self.reactions_bounds.append(ExactSubstrateReactionBounds(substrate=d_glucose,bound=-1))
        if target==ATPM:
            self.reactions_bounds.append(ReactionBounds(reaction_name=ATPM,lower_bound=0))
            self.objective_function_name=ATPM
            self.target=None
        else:
            self.reactions_bounds.append(ReactionBounds(reaction_name=ATPM,lower_bound=0,upper_bound=0))
            self.objective_function_name=target+CofactorAndPrecursorsFBATest.drain_suffix
            self.target=target.lower()
            
    def run(self,model):
        if self.target:
            oxidized_target=get_cytosolic_metabolite(model,self.target[:-1])
            reduced_target=get_cytosolic_metabolite(model,self.target)
            released_proton=get_cytosolic_metabolite(model,proton)
            drain_reaction=cobra.Reaction(self.target.upper()+CofactorAndPrecursorsFBATest.drain_suffix)
            drain_reaction.name=self.target+" drain reaction"
            drain_reaction.add_metabolites({reduced_target:-1.0,
                                            oxidized_target:1.0,
                                            released_proton:1.0})
            model.add_reaction(drain_reaction)
        opt_solution=super(CofactorAndPrecursorsFBATest,self).run(model)
        return opt_solution
    
    def run_experiment_and_print_result(self,model):
        opt_solution=self.run(model)
        print(self.name+": "+
            str(round(opt_solution.objective_value,4))+"mol/mol_glucose"+
            "  ATP shadow price: "+str(opt_solution.shadow_prices[ADP+cytosolic_metabolite_suffix]))
        
def get_cytosolic_metabolite(model,metabolite):
    return model.metabolites.get_by_id(metabolite+cytosolic_metabolite_suffix)
    
class ReactionBounds:
    def __init__(self,reaction_name,lower_bound=None,upper_bound=None):
        self.reaction_name=reaction_name
        self.upper_bound=upper_bound
        self.lower_bound=lower_bound
        
    def set_reaction_bounds(self,model):
        reaction=model.reactions.get_by_id(self.reaction_name)
        if not self.lower_bound==None:
            reaction.lower_bound=self.lower_bound
        if not self.upper_bound==None:
            reaction.upper_bound=self.upper_bound

class RealisticSubstrateReactionBounds(ReactionBounds):
    def __init__(self,substrate):
        super(RealisticSubstrateReactionBounds,self).__init__(exchange_reaction_name_for_substrate(substrate),
                                                                lower_bound=organic_exchange_realistic_lower_bound if not substrate==d_glucose else glucose_exchange_realistic_lower_bound,
                                                                upper_bound=unlimited_reaction_bound)

class UnlimitedSubstrateReactionBounds(ReactionBounds):
    def __init__(self,substrate):
        super(UnlimitedSubstrateReactionBounds,self).__init__(exchange_reaction_name_for_substrate(substrate),
                                                                lower_bound=-unlimited_reaction_bound,
                                                                upper_bound=unlimited_reaction_bound)

class UnavailableSubstrateReactionBounds(ReactionBounds):
    def __init__(self,substrate):
        super(UnavailableSubstrateReactionBounds,self).__init__(exchange_reaction_name_for_substrate(substrate),
                                                                lower_bound=limited_reaction_bound,
                                                                upper_bound=unlimited_reaction_bound)
                                                                
class ExactSubstrateReactionBounds(ReactionBounds):
    def __init__(self,substrate,bound):
        super(ExactSubstrateReactionBounds,self).__init__(exchange_reaction_name_for_substrate(substrate),
                                                                lower_bound=bound,
                                                                upper_bound=bound)
    
def run_experiments_and_print_result(model,experiments):
    for experiment in experiments:
        with model as model:
            experiment.run_experiment_and_print_result(model)

def exchange_reaction_name_for_substrate(substrateName):
    return "EX_"+substrateName+"_e"

def run_example_1():
    print("Maximal growth rate on organic substrates:")
    model=cobra.test.create_test_model("textbook")
    substrates=[acetate,acetaldehyde,two_oxoglutarate,ethanol,d_fructose,fumarate,d_glucose,l_glutamine,l_glutamate,d_lactate,l_malate,pyruvate,succinate]
    experiments=[]
    for substrate in substrates:
        for aerobic in [True,False]:
            experiments.append(SubstrateFBAExperiment(substrate,aerobic))
    run_experiments_and_print_result(model,experiments)
    
def run_example_2():
    print("Maximum yield of cofactors and precursors")
    model=cobra.test.create_test_model("textbook")
    experiments=[]
    for cofactor in [ATPM,"NADH","NADPH"]:
        for aerobic in [True,False]:
            experiments.append(CofactorAndPrecursorsFBATest(cofactor,aerobic))
    run_experiments_and_print_result(model,experiments)
    
def run_example_3():
    model=cobra.test.create_test_model("textbook")
    target=get_cytosolic_metabolite(model,"3pg")
    drain_reaction=cobra.Reaction(target.id.upper()+CofactorAndPrecursorsFBATest.drain_suffix)
    drain_reaction.name=target.id+" drain reaction"
    drain_reaction.add_metabolites({target:-1.0})
    model.add_reaction(drain_reaction)
    model.objective=drain_reaction
    print(model.optimize().objective_value)
    
def run_example_4():
    model=cobra.test.create_test_model("textbook")
    model.reactions.get_by_id("EX_glc__D_e").lower_bound=0
    model.reactions.get_by_id("EX_succ_e").lower_bound=-20
    model.reactions.get_by_id("EX_succ_e").upper_bound=-20
    solution=model.optimize()
    model.reactions.get_by_id("Biomass_Ecoli_core").lower_bound=solution.f
    model.reactions.get_by_id("Biomass_Ecoli_core").upper_bound=solution.f
    model.objective=model.reactions.get_by_id("ME1")
    min_solution=model.optimize(objective_sense="minimize")
    max_solution=model.optimize(objective_sense="maximize")
    print(min_solution.objective_value)
    print(max_solution.objective_value)
    
def run_example_5():
    model=cobra.test.create_test_model("textbook")
    print(model.fluxVariability)