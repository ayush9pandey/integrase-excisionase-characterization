from biocrnpyler import *

from itertools import chain
class IntegraseSite(Component):
    def __init__(self, integrase: Species, 
                 curr_dna: Species, reversed_dna: Species, 
                 name: str=None, cooperativity: int=1, parameter_dict = None,
                 mechanism = None, **keywords):
        if name is None:
            name = integrase.name
        self.name = name
        self.integrase = integrase
        self.curr_dna = curr_dna
        self.reversed_dna = reversed_dna
        self.complex_species = None
        self.flipped_complex = None
        self.species = None
        self.cooperativity = cooperativity
        self.parameter_dict = parameter_dict
        mechanisms = {mechanism.mechanism_type:mechanism}
        Component.__init__(self=self, name=name, mechanisms=mechanisms, **keywords)
        
    def get_species(self) -> List[Species]:
        return self.species

    def update_species(self) -> List[Species]:
        integrase_flipper = self.get_mechanism('integrase-flipper')
        self.species = integrase_flipper.update_species(integrase = self.integrase, 
                                                        curr_dna = self.curr_dna, 
                                                        reversed_dna = self.reversed_dna, 
                                                        cooperativity = self.cooperativity)
        return self.species

    def update_reactions(self) -> List[Reaction]:
        integrase_flipper = self.get_mechanism('integrase-flipper')
        return integrase_flipper.update_reactions(integrase = self.integrase, curr_dna = self.curr_dna, 
                                                  reversed_dna = self.reversed_dna, cooperativity = self.cooperativity,
                                                  complex_species = self.complex_species, 
                                                  flipped_complex = self.flipped_complex, 
                                                  component = self, 
                                                  part_id = None, parameter_dict = self.parameter_dict)


class IntegraseFlipper(Mechanism):
    """Mechanism for the simple integrase flipping. 
    Integrase binds to attP-attB site, forms a complex,
    the complex flips the site to attL-attR, then integrase unbinds.
    """
    def __init__(self, name: str="integrase-flipper", mechanism_type: str="integrase-flipper", **keywords):
        """Initializes a IntegraseFlipper instance.
        :param name: name of the Mechanism, default: integrase-flipper
        :param mechanism_type: type of the Mechanism, default: integrase-flipper
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)

    def update_species(self, integrase, curr_dna, reversed_dna, cooperativity):
        list_of_species = [cooperativity*[integrase], [curr_dna]]
        # attP_attB_4I
        complex_species = Complex(species = list_of_species)
        self.complex_species = complex_species
        # attL_attR_4I
        list_of_species = [cooperativity*[integrase], [reversed_dna]]
        flipped_complex = Complex(species = list_of_species)
        self.flipped_complex = flipped_complex
        return [integrase, curr_dna, reversed_dna, 
                complex_species, flipped_complex]

    def update_reactions(self, integrase, curr_dna, reversed_dna, cooperativity,
                         complex_species, flipped_complex, component = None, 
                         part_id = None, parameter_dict = None):
        if part_id is None and component is not None:
            part_id = component.name
        if parameter_dict is None and component is None:
            raise ValueError("Must pass in either a component or a parameter dictionary.")
        elif parameter_dict is None:
            a_i = component.get_parameter("a_i", part_id = part_id, mechanism = self)
            d_i = component.get_parameter("d_i", part_id = part_id, mechanism = self)
            k_int = component.get_parameter("k_int", part_id = part_id, mechanism = self)
            a_r = component.get_parameter("a_r", part_id = part_id, mechanism = self)
            d_r = component.get_parameter("d_r", part_id = part_id, mechanism = self)
        else:
            a_i = parameter_dict['a_i']
            d_i = parameter_dict['d_i']
            k_int = parameter_dict['k_int']
            a_r = parameter_dict['a_r']
            d_r = parameter_dict['d_r']
        
        if complex_species is None:
            if self.complex_species is None:
                raise ValueError('The complex_species is not set for mechanism {0}'.format(self.name))
            complex_species = self.complex_species
        if flipped_complex is None:
            if self.flipped_complex is None:
                raise ValueError('The flipped_species is not set for mechanism {0}'.format(self.name))
            flipped_complex = self.flipped_complex   
        list_of_species = [cooperativity*[integrase], [curr_dna]]
        integrase_binding_reaction = Reaction.from_massaction([item for sublist in list_of_species for item in sublist],
                                                              [complex_species],
                                                              k_forward = a_i, k_reverse = d_i)
        flipping_reaction = Reaction.from_massaction([complex_species], [flipped_complex], k_forward = k_int)
        
        list_of_species = [cooperativity*[integrase], [reversed_dna]]
        integrase_unbinding_reaction = Reaction.from_massaction([flipped_complex],
                                                                [item for sublist in list_of_species for item in sublist],
                                                                k_forward = a_r, k_reverse = d_r)
        reactions = [integrase_binding_reaction, flipping_reaction, integrase_unbinding_reaction]
        return reactions

from itertools import chain
class ExcisionaseSite(Component):
    def __init__(self, excisionase: Species, integrase_site: Component,
                 curr_dna: Species, reversed_dna: Species,
                 name: str=None, cooperativity: int=1, parameter_dict = None,
                 mechanism = None, **keywords):

        if name is None:
            name = excisionase.name
        self.name = name
        self.excisionase = excisionase # E
        self.curr_dna = curr_dna # A 
        self.reversed_dna = reversed_dna # R
        integrase_site.update_species()
        self.integrase_site = integrase_site
        self.complex_species = None
        self.flipped_complex = None
        self.integrase_excisionase_complex = None
        self.species = None
        self.cooperativity = cooperativity
        self.parameter_dict = parameter_dict
        mechanisms = {mechanism.mechanism_type:mechanism}
        self.mechanism_type = mechanism.mechanism_type
        Component.__init__(self=self, name=name, mechanisms=mechanisms, **keywords)
        
    def get_species(self) -> List[Species]:
        return self.species

    def update_species(self) -> List[Species]:
        excisionase_flipper = self.get_mechanism(self.mechanism_type)
        self.species = excisionase_flipper.update_species(excisionase = self.excisionase, 
                                                          curr_dna = self.curr_dna, 
                                                          reversed_dna = self.reversed_dna, 
                                                          cooperativity = self.cooperativity,
                                                          integrase_site = self.integrase_site)
        return self.species
        

    def update_reactions(self) -> List[Reaction]:
        xis_flipper = self.get_mechanism(self.mechanism_type)
        return xis_flipper.update_reactions(excisionase = self.excisionase, 
                                            curr_dna = self.curr_dna, 
                                            reversed_dna = self.reversed_dna, 
                                            cooperativity = self.cooperativity,
                                            integrase_site = self.integrase_site, 
                                            complex_species = self.complex_species,
                                            flipped_complex = self.flipped_complex, 
                                            integrase_excisionase_complex = self.integrase_excisionase_complex,
                                            component = self, part_id = None, parameter_dict = self.parameter_dict)
   
    
class ExcisionaseFlipper(Mechanism):
    """Mechanism for the excisionase flipping. 
    Excisionase binds to attL-attR site, forms a complex,
    the complex flips the site to attP-attB, then excisionase unbinds.
    Both binding of excisionase directly to the integrase-dna complex and
    first binding to integrase then to DNA are allowed.
    """
    def __init__(self, name: str="excisionase-flipper", 
                 mechanism_type: str="excisionase-flipper", **keywords):
        """Initializes a ExcisionaseFlipper instance.
        :param name: name of the Mechanism, default: excisionase-flipper
        :param mechanism_type: type of the Mechanism, default: excisionase-flipper
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)   
    
    def update_species(self, excisionase, curr_dna, reversed_dna, cooperativity, integrase_site):
        # Update integrase_site species
        i1, i2, i3, i4, i5 = integrase_site.update_species()
        integrase_site.integrase = i1
        integrase_site.curr_dna = i2
        integrase_site.reversed_dna = i3
        integrase_site.complex_species = i4
        integrase_site.flipped_complex = i5
        # Excisionase binds to the DNA with integrase already bound to it...
        list_of_species = [cooperativity*[excisionase], 
                           [integrase_site.flipped_complex]]
#         attL_attR_4I_4E
        complex_species = Complex(species = list_of_species) # C_E
        self.complex_species = complex_species
#         attP_attB_4I_4E
        list_of_species = [cooperativity*[excisionase], [integrase_site.complex_species]]
        flipped_complex = Complex(list_of_species) # C_R
        self.flipped_complex = flipped_complex
        # ... or it binds to integrase first then binds the DNA.
        integrase_excisionase_complex = Complex(species = [integrase_site.integrase, excisionase]) # C_IE
        self.integrase_excisionase_complex = integrase_excisionase_complex
        
        return [excisionase, curr_dna, reversed_dna, complex_species, 
                flipped_complex, integrase_excisionase_complex]

    def update_reactions(self, excisionase, curr_dna, reversed_dna, cooperativity,
                         integrase_site, complex_species, flipped_complex, integrase_excisionase_complex,
                         component = None, part_id = None, parameter_dict = None):
        if part_id is None and component is not None:
            part_id = component.name
        if parameter_dict is None and component is None:
            raise ValueError("Must pass in either a component or a parameter dictionary.")
        elif parameter_dict is None:
            a_e1 = component.get_parameter("a_e1", part_id = part_id, mechanism = self)
            d_e1 = component.get_parameter("d_e1", part_id = part_id, mechanism = self)
            a_r1 = component.get_parameter("a_r1", part_id = part_id, mechanism = self)
            d_r1 = component.get_parameter("d_r1", part_id = part_id, mechanism = self)
            a_s1 = component.get_parameter("a_s1", part_id = part_id, mechanism = self)
            d_s1 = component.get_parameter("d_s1", part_id = part_id, mechanism = self)
            a_s2 = component.get_parameter("a_s2", part_id = part_id, mechanism = self)
            d_s2 = component.get_parameter("d_s2", part_id = part_id, mechanism = self)
            a_e2 = component.get_parameter("a_e2", part_id = part_id, mechanism = self)
            d_e2 = component.get_parameter("d_e2", part_id = part_id, mechanism = self)
            a_r2 = component.get_parameter("a_r2", part_id = part_id, mechanism = self)
            d_r2 = component.get_parameter("d_r2", part_id = part_id, mechanism = self)
            k_exc = component.get_parameter("k_exc", part_id = part_id, mechanism = self)
        else:
            a_e1 = parameter_dict["a_e1"]
            d_e1 = parameter_dict["d_e1"]
            a_r1 = parameter_dict["a_r1"]
            d_r1 = parameter_dict["d_r1"]
            a_s1 = parameter_dict["a_s1"]
            d_s1 = parameter_dict["d_s1"]
            a_s2 = parameter_dict["a_s2"]
            d_s2 = parameter_dict["d_s2"]
            a_e2 = parameter_dict["a_e2"]
            d_e2 = parameter_dict["d_e2"]
            a_r2 = parameter_dict["a_r2"]
            d_r2 = parameter_dict["d_r2"]
            k_exc = parameter_dict["k_exc"]
        if complex_species is None:
            if self.complex_species is None:
                raise ValueError('The complex_species is not set for mechanism {0}'.format(self.name))
            complex_species = self.complex_species
        if flipped_complex is None:
            if self.flipped_complex is None:
                raise ValueError('The flipped_species is not set for mechanism {0}'.format(self.name))
            flipped_complex = self.flipped_complex
        if integrase_excisionase_complex is None:
            if self.integrase_excisionase_complex is None:
                raise ValueError('The integrase_excisionase_species is not set for mechanism {0}'.format(self.name))
            integrase_excisionase_complex = self.integrase_excisionase_complex
        # First set of excisionase action
        list_of_species = [cooperativity*[excisionase], [integrase_site.flipped_complex]]
        excisionase_binding_reaction = Reaction.from_massaction([item for sublist in list_of_species 
                                                                 for item in sublist],
                                                                [complex_species],
                                                                k_forward = a_e1, k_reverse = d_e1)
        flipping_reaction = Reaction.from_massaction([complex_species], 
                                                     [flipped_complex], k_forward = k_exc)
        
        list_of_species = [cooperativity*[excisionase],
                           integrase_site.cooperativity*[integrase_site.integrase], 
                           [reversed_dna]]
        excisionase_unbinding_reaction = Reaction.from_massaction([flipped_complex],
                                                                [item for sublist in list_of_species 
                                                                 for item in sublist],
                                                                k_forward = a_r1, k_reverse = d_r1)
        # Excisionase sequestration reactions
        integrase_excisionase_bind = Reaction.from_massaction([integrase_site.integrase, excisionase],
                                                             [integrase_excisionase_complex],
                                                             k_forward = a_s1, k_reverse = d_s1)
        
        list_of_species = [cooperativity*[excisionase], [integrase_site.complex_species]]
        excisionase_binding_reaction2 = Reaction.from_massaction([item for sublist in list_of_species 
                                                                 for item in sublist],
                                                                 [flipped_complex],
                                                                 k_forward = a_s2, k_reverse = d_s2)
        # Second set of excisionase action
        list_of_species = [cooperativity*[integrase_excisionase_complex], [curr_dna]]
        int_exc_binding_curr_dna = Reaction.from_massaction([item for sublist in list_of_species 
                                                                 for item in sublist],
                                                                 [complex_species],
                                                                 k_forward = a_e2, k_reverse = d_e2)
        list_of_species = [cooperativity*[integrase_excisionase_complex], [reversed_dna]]
        int_exc_binding_rev_dna = Reaction.from_massaction([item for sublist in list_of_species 
                                                                 for item in sublist],
                                                                 [flipped_complex],
                                                                 k_forward = a_r2, k_reverse = d_r2)
        
        reactions = [excisionase_binding_reaction, flipping_reaction, excisionase_unbinding_reaction,
                     excisionase_binding_reaction2, integrase_excisionase_bind, int_exc_binding_curr_dna,
                     int_exc_binding_rev_dna]
        return reactions

class SimpleExcisionaseFlipper(Mechanism):
    """Mechanism for the simple excisionase flipping. 
    Excisionase binds to integrase then to attL-attR site, forms a complex,
    the complex flips the site to attP-attB, then excisionase and integrase unbind.
    Other kinds of interactions are not allowed.
    """
    def __init__(self, name: str="simple-excisionase-flipper", 
                 mechanism_type: str="simple-excisionase-flipper", **keywords):
        """Initializes a SimpleExcisionaseFlipper instance.
        :param name: name of the Mechanism, default: simple-excisionase-flipper
        :param mechanism_type: type of the Mechanism, default: simple-excisionase-flipper
        :param keywords:
        """
        Mechanism.__init__(self, name, mechanism_type)   
    
    def update_species(self, excisionase, curr_dna, reversed_dna, cooperativity, integrase_site):
        # Update integrase_site species
        i1, i2, i3, i4, i5 = integrase_site.update_species()
        integrase_site.integrase = i1
        integrase_site.curr_dna = i2
        integrase_site.reversed_dna = i3
        integrase_site.complex_species = i4
        integrase_site.flipped_complex = i5
        # Excisionase binds to the DNA with integrase already bound to it...
        list_of_species = [cooperativity*[excisionase], 
                           [integrase_site.flipped_complex]]
#         attL_attR_4I_4E
        complex_species = Complex(species = list_of_species) # C_E
        self.complex_species = complex_species
#         attP_attB_4I_4E
        list_of_species = [cooperativity*[excisionase], [integrase_site.complex_species]]
        flipped_complex = Complex(list_of_species) # C_R
        self.flipped_complex = flipped_complex
        # ... or it binds to integrase first then binds the DNA.
        integrase_excisionase_complex = Complex(species = [integrase_site.integrase, excisionase]) # C_IE
        self.integrase_excisionase_complex = integrase_excisionase_complex
        
        return [excisionase, curr_dna, reversed_dna, complex_species, 
                flipped_complex, integrase_excisionase_complex]

    def update_reactions(self, excisionase, curr_dna, reversed_dna, cooperativity,
                         integrase_site, complex_species, flipped_complex, integrase_excisionase_complex,
                         component = None, part_id = None, parameter_dict = None):
        if part_id is None and component is not None:
            part_id = component.name
        if parameter_dict is None and component is None:
            raise ValueError("Must pass in either a component or a parameter dictionary.")
        elif parameter_dict is None:
            a_e1 = component.get_parameter("a_e1", part_id = part_id, mechanism = self)
            d_e1 = component.get_parameter("d_e1", part_id = part_id, mechanism = self)
            a_r1 = component.get_parameter("a_r1", part_id = part_id, mechanism = self)
            d_r1 = component.get_parameter("d_r1", part_id = part_id, mechanism = self)
            a_s1 = component.get_parameter("a_s1", part_id = part_id, mechanism = self)
            d_s1 = component.get_parameter("d_s1", part_id = part_id, mechanism = self)
            a_s2 = component.get_parameter("a_s2", part_id = part_id, mechanism = self)
            d_s2 = component.get_parameter("d_s2", part_id = part_id, mechanism = self)
            a_e2 = component.get_parameter("a_e2", part_id = part_id, mechanism = self)
            d_e2 = component.get_parameter("d_e2", part_id = part_id, mechanism = self)
            a_r2 = component.get_parameter("a_r2", part_id = part_id, mechanism = self)
            d_r2 = component.get_parameter("d_r2", part_id = part_id, mechanism = self)
            k_exc = component.get_parameter("k_exc", part_id = part_id, mechanism = self)
        else:
            a_e1 = parameter_dict["a_e1"]
            d_e1 = parameter_dict["d_e1"]
            a_r1 = parameter_dict["a_r1"]
            d_r1 = parameter_dict["d_r1"]
            a_s1 = parameter_dict["a_s1"]
            d_s1 = parameter_dict["d_s1"]
            a_s2 = parameter_dict["a_s2"]
            d_s2 = parameter_dict["d_s2"]
            a_e2 = parameter_dict["a_e2"]
            d_e2 = parameter_dict["d_e2"]
            a_r2 = parameter_dict["a_r2"]
            d_r2 = parameter_dict["d_r2"]
            k_exc = parameter_dict["k_exc"]
        if complex_species is None:
            if self.complex_species is None:
                raise ValueError('The complex_species is not set for mechanism {0}'.format(self.name))
            complex_species = self.complex_species
        if flipped_complex is None:
            if self.flipped_complex is None:
                raise ValueError('The flipped_species is not set for mechanism {0}'.format(self.name))
            flipped_complex = self.flipped_complex
        if integrase_excisionase_complex is None:
            if self.integrase_excisionase_complex is None:
                raise ValueError('The integrase_excisionase_species is not set for mechanism {0}'.format(self.name))
            integrase_excisionase_complex = self.integrase_excisionase_complex
        # First set of excisionase action
#         list_of_species = [cooperativity*[excisionase], [integrase_site.flipped_complex]]
#         excisionase_binding_reaction = Reaction.from_massaction([item for sublist in list_of_species 
#                                                                  for item in sublist],
#                                                                 [complex_species],
#                                                                 k_forward = a_e1, k_reverse = d_e1)
        flipping_reaction = Reaction.from_massaction([complex_species], 
                                                     [flipped_complex], k_forward = k_exc)
        
        list_of_species = [cooperativity*[excisionase],
                           integrase_site.cooperativity*[integrase_site.integrase], 
                           [reversed_dna]]
        excisionase_unbinding_reaction = Reaction.from_massaction([flipped_complex],
                                                                [item for sublist in list_of_species 
                                                                 for item in sublist],
                                                                k_forward = a_r1, k_reverse = d_r1)
        # Excisionase sequestration reactions
        integrase_excisionase_bind = Reaction.from_massaction([integrase_site.integrase, excisionase],
                                                             [integrase_excisionase_complex],
                                                             k_forward = a_s1, k_reverse = d_s1)
        
#         list_of_species = [cooperativity*[excisionase], [integrase_site.complex_species]]
#         excisionase_binding_reaction2 = Reaction.from_massaction([item for sublist in list_of_species 
#                                                                  for item in sublist],
#                                                                  [flipped_complex],
#                                                                  k_forward = a_s2, k_reverse = d_s2)
        # Second set of excisionase action
        list_of_species = [cooperativity*[integrase_excisionase_complex], [curr_dna]]
        int_exc_binding_curr_dna = Reaction.from_massaction([item for sublist in list_of_species 
                                                                 for item in sublist],
                                                                 [complex_species],
                                                                 k_forward = a_e2, k_reverse = d_e2)
        list_of_species = [cooperativity*[integrase_excisionase_complex], [reversed_dna]]
        int_exc_binding_rev_dna = Reaction.from_massaction([item for sublist in list_of_species 
                                                                 for item in sublist],
                                                                 [flipped_complex],
                                                                 k_forward = a_r2, k_reverse = d_r2)
        
        reactions = [flipping_reaction, excisionase_unbinding_reaction,
                     integrase_excisionase_bind, int_exc_binding_curr_dna,
                     int_exc_binding_rev_dna]
        return reactions

class TxTlMixture(Mixture):
    def __init__(self, name="", rnap="RNAP", ribosome="Ribo", rnaase="RNAase", **kwargs):
        Mixture.__init__(self, name=name, **kwargs)

        # Create Components for TxTl machinery
        self.rnap = Protein(rnap)
        self.ribosome = Protein(ribosome)
        self.rnaase = Protein(rnaase)
        default_components = [
            self.rnap, self.ribosome, self.rnaase
        ]
        self.add_components(default_components)

        #Create TxTl Mechansisms
        mech_tx = Transcription_MM(rnap = self.rnap.get_species())
        mech_tl = Translation_MM(ribosome = self.ribosome.get_species())
        mech_rna_deg = Degredation_mRNA_MM(nuclease=self.rnaase.get_species())
        mech_cat = MichaelisMenten()
        mech_bind = One_Step_Binding()

       #Create Global Dilution Mechanisms
        dilution_mechanism = Dilution(filter_dict = {"rna":True, "machinery":False}, 
                                      default_on = False)
        default_mechanisms = {
            mech_tx.mechanism_type: mech_tx,
            mech_tl.mechanism_type: mech_tl,
            mech_cat.mechanism_type: mech_cat,
            mech_bind.mechanism_type:mech_bind,
            "dilution":mech_rna_deg,
        }

        self.add_mechanisms(default_mechanisms)