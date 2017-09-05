import cobra.test
import warnings
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from cobra.flux_analysis import flux_variability_analysis

from matplotlib.colors import LightSource
from matplotlib import cm


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
malic_enzyme_NAD="ME1"
malic_enzyme_NADPH="ME2"
fumarate_reductaase="FRD7"
succinate_dehydrogenase="SUCDi"
malate_dehydrogenase="MDH"
NAD_transdehydrogenase="NADTRHD"
pep_carboxykinase="PPCK"
pyruvate_kinase="PYK"

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
        if aerobic==None:
            pass
        elif aerobic==True:
            self.name=self.name+"_aerobic"
            self.reactions_bounds.append(UnlimitedSubstrateReactionBounds(oxygen))
        elif aerobic==False:
            self.name=self.name+"_unaerobic"
            self.reactions_bounds.append(UnavailableSubstrateReactionBounds(oxygen)) 
        
    def run(self,model):
        return self.set_model_bounds_and_get_optimal_values(model)
        
    def set_model_bounds(self,model):
        for reaction_bound in self.reactions_bounds:
            reaction_bound.set_reaction_bounds(model)
        
    def set_model_bounds_and_get_optimal_values(self,model):
        self.set_model_bounds(model)
        if self.objective_function_name:
            model.objective=model.reactions.get_by_id(self.objective_function_name)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            solution=model.optimize("max")
            return(solution)
      
    def run_experiment_and_print_result(self,model):
        pass
    
            
class SubstrateFBAExperiment(FBAExperiment):
    def __init__(self,substrate,aerobic,specific_substrate_bound=None):
        super(SubstrateFBAExperiment,self).__init__(name="Realistic_"+substrate,aerobic=aerobic)
        if specific_substrate_bound==None:
            self.reactions_bounds.append(RealisticSubstrateReactionBounds(substrate))
        else:
            self.reactions_bounds.append(ExactSubstrateReactionBounds(substrate,specific_substrate_bound))
        if not substrate==d_glucose:
            self.reactions_bounds.append(UnavailableSubstrateReactionBounds(d_glucose))
    
    def run(self,model):
        return super(SubstrateFBAExperiment,self).run(model)
        
    def run_experiment_and_print_result(self,model):
        opt_solution=self.run(model)
        opt_value=opt_solution.objective_value if opt_solution.status=="optimal" else 0.0
        print(self.name+": "+str(round(opt_value,4))+"/h")
        
class SubstrateAndOxygenLimitedFBAExperiment(SubstrateFBAExperiment):
    def __init__(self,substrate,exact_substrate_bound,exact_oxygen_bound):
        super(SubstrateAndOxygenLimitedFBAExperiment,self).__init__(substrate=substrate,aerobic=None,specific_substrate_bound=exact_substrate_bound)
        self.reactions_bounds.append(ExactSubstrateReactionBounds(oxygen,exact_oxygen_bound))
        
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
            
class FluxVariabilityExperiment(SubstrateFBAExperiment):
    def __init__(self,substrate,reactions_to_test,aerobic,specific_substrate_bound):
        super(FluxVariabilityExperiment,self).__init__(substrate,aerobic,specific_substrate_bound)
        self.reactions_to_test=list(reactions_to_test)
    
    def run(self,model):
        super(FluxVariabilityExperiment,self).set_model_bounds(model)
        return(flux_variability_analysis(model, self.reactions_to_test))
        
    def run_experiment_and_print_result(self,model):
        fva_result=self.run(model)
        print(fva_result)
        
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
    
def run_example_2_1():
    print("Maximum yield of cofactors and precursors")
    model=cobra.test.create_test_model("textbook")
    target=get_cytosolic_metabolite(model,"3pg")
    drain_reaction=cobra.Reaction(target.id.upper()+CofactorAndPrecursorsFBATest.drain_suffix)
    drain_reaction.name=target.id+" drain reaction"
    drain_reaction.add_metabolites({target:-1.0})
    model.add_reaction(drain_reaction)
    model.objective=drain_reaction
    print(model.optimize().objective_value)
    
def run_example_3_1():
    print("Alternate optimal solutions")
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
    
def run_example_3():
    print("Alternate optimal solutions")
    model=cobra.test.create_test_model("textbook")
    reactions_to_test=[model.reactions.get_by_id(reaction) for reaction in {malic_enzyme_NAD,
                                                                            malic_enzyme_NADPH,
                                                                            fumarate_reductaase,
                                                                            succinate_dehydrogenase,
                                                                            malate_dehydrogenase,
                                                                            NAD_transdehydrogenase,
                                                                            pep_carboxykinase,
                                                                            pyruvate_kinase}]
    """
    model.reactions.get_by_id("EX_glc__D_e").lower_bound=0
    model.reactions.get_by_id("EX_succ_e").lower_bound=-20
    model.reactions.get_by_id("EX_succ_e").upper_bound=-20
    reactions_to_test=[model.reactions.get_by_id(reaction) for reaction in {malic_enzyme_NAD,malic_enzyme_NADPH,fumarate_reductaase,succinate_dehydrogenase,malate_dehydrogenase,NAD_transdehydrogenase,pep_carboxykinase,pyruvate_kinase}]
    print(flux_variability_analysis(model, reactions_to_test))
    """
    experiment=FluxVariabilityExperiment(substrate=succinate,
                                        reactions_to_test=reactions_to_test,
                                        aerobic=True,
                                        specific_substrate_bound=-20.0)
    experiment.run_experiment_and_print_result(model)
    
def run_example_4_1():
    print("Robustness analysis, variable glucose")
    model=cobra.test.create_test_model("textbook")
    """
    oxygen_exchange_reaction=model.reactions.get_by_id(exchange_reaction_name_for_substrate(oxygen))
    oxygen_exchange_reaction.lower_bound=-17
    oxygen_exchange_reaction.upper_bound=-17
    growth_rates=[]
    for i in range(0,21):
        with model as model:
            glc_ex_reaction=model.reactions.get_by_id(exchange_reaction_name_for_substrate(d_glucose))
            glc_ex_reaction.lower_bound=-i
            glc_ex_reaction.upper_bound=-i
            growth_rates.append(model.optimize().objective_value)
    plt.plot(growth_rates)
    plt.show()
    """
    growth_rates=[]
    for i in range(0,21):
        with model as model:
            growth_rates.append(SubstrateAndOxygenLimitedFBAExperiment(substrate=d_glucose,exact_substrate_bound=-i,exact_oxygen_bound=-17).run(model).objective_value)
    plt.plot(growth_rates)
    plt.show()
    
def run_example_4_2():
    print("Robustness analysis, variable oxygen")
    model=cobra.test.create_test_model("textbook")
    """
    glc_ex_reaction=model.reactions.get_by_id(exchange_reaction_name_for_substrate(d_glucose))
    glc_ex_reaction.lower_bound=-10
    glc_ex_reaction.upper_bound=-10
    growth_rates=[]
    for i in range(0,26):
        with model as model:
            oxygen_ex_reaction=model.reactions.get_by_id(exchange_reaction_name_for_substrate(oxygen))
            oxygen_ex_reaction.lower_bound=-i
            oxygen_ex_reaction.upper_bound=-i
            growth_rates.append(model.optimize().objective_value)
    plt.plot(growth_rates)
    plt.show()
    """
    growth_rates=[]
    for i in range(0,26):
        with model as model:
            growth_rates.append(SubstrateAndOxygenLimitedFBAExperiment(substrate=d_glucose,exact_substrate_bound=-10,exact_oxygen_bound=-i).run(model).objective_value)
    plt.plot(growth_rates)
    plt.show()
    
def run_example_5():
    print("Phenotypic phase analysis")
    model=cobra.test.create_test_model("textbook")
    glucose_upper_limit=21
    oxygen_upper_limit=21
    glucoses = np.arange(glucose_upper_limit)
    oxygens = np.arange(oxygen_upper_limit)
    glucoses, oxygens = np.meshgrid(glucoses, oxygens)
    growth_rates=calc_growth_rate(glucoses,oxygens,model)
    
    plt.clf()
    
    light = LightSource(45, 60)
    
    white = np.ones((growth_rates.shape[0], growth_rates.shape[1], 3))
    green = white * np.array([0,1,0])
    
    
    illuminated_surface = light.shade_rgb(green,growth_rates)
    ax = plt.axes(projection='3d')
    ax.plot_surface(glucoses, oxygens, growth_rates, cstride=1, rstride=1, facecolors=illuminated_surface)
    plt.show()
    plt.clf()
    
    plt.imshow(growth_rates, cmap='hot', interpolation='nearest')
    plt.show()
    plt.clf()
    
    
@np.vectorize
def calc_growth_rate(glucose_bound,oxygen_bound,model):
    with model as model:
        fba_solution=SubstrateAndOxygenLimitedFBAExperiment(substrate=d_glucose,exact_substrate_bound=-glucose_bound,exact_oxygen_bound=-oxygen_bound).run(model)
        return fba_solution.objective_value if fba_solution.status=="optimal" else 0.0