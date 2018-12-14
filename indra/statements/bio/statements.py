"""
Statements represent mechanistic relationships between biological agents.

Statement classes follow an inheritance hierarchy, with all Statement types
inheriting from the parent class :py:class:`Statement`. At
the next level in the hierarchy are the following classes:

- :py:class:`Complex`
- :py:class:`Modification`
- :py:class:`SelfModification`
- :py:class:`RegulateActivity`
- :py:class:`RegulateAmount`
- :py:class:`ActiveForm`
- :py:class:`Translocation`
- :py:class:`Gef`
- :py:class:`Gap`
- :py:class:`Conversion`

There are several types of Statements representing post-translational
modifications that further inherit from
:py:class:`Modification`:

- :py:class:`Phosphorylation`
- :py:class:`Dephosphorylation`
- :py:class:`Ubiquitination`
- :py:class:`Deubiquitination`
- :py:class:`Sumoylation`
- :py:class:`Desumoylation`
- :py:class:`Hydroxylation`
- :py:class:`Dehydroxylation`
- :py:class:`Acetylation`
- :py:class:`Deacetylation`
- :py:class:`Glycosylation`
- :py:class:`Deglycosylation`
- :py:class:`Farnesylation`
- :py:class:`Defarnesylation`
- :py:class:`Geranylgeranylation`
- :py:class:`Degeranylgeranylation`
- :py:class:`Palmitoylation`
- :py:class:`Depalmitoylation`
- :py:class:`Myristoylation`
- :py:class:`Demyristoylation`
- :py:class:`Ribosylation`
- :py:class:`Deribosylation`
- :py:class:`Methylation`
- :py:class:`Demethylation`

There are additional subtypes of :py:class:`SelfModification`:

- :py:class:`Autophosphorylation`
- :py:class:`Transphosphorylation`

Interactions between proteins are often described simply in terms of their
effect on a protein's "activity", e.g., "Active MEK activates ERK", or "DUSP6
inactives ERK".  These types of relationships are indicated by the
:py:class:`RegulateActivity` abstract base class which has subtypes

- :py:class:`Activation`
- :py:class:`Inhibition`

while the :py:class:`RegulateAmount` abstract base class has subtypes

- :py:class:`IncreaseAmount`
- :py:class:`DecreaseAmount`

Statements involve one or more *Concepts*, which, depending on the
semantics of the Statement, are typically biological *Agents*,
such as proteins, represented by the class :py:class:`Agent`.
Agents can have several types of context specified on them including

- a specific post-translational modification state (indicated by one or
  more instances of :py:class:`ModCondition`),
- other bound Agents (:py:class:`BoundCondition`),
- mutations (:py:class:`MutCondition`),
- an activity state (:py:class:`ActivityCondition`), and
- cellular location

The *active* form of an agent (in terms of its post-translational modifications
or bound state) is indicated by an instance of the class
:py:class:`ActiveForm`.

Agents also carry grounding information which links them to database entries.
These database references are represented as a dictionary in the `db_refs`
attribute of each Agent. The dictionary can have multiple entries. For
instance, INDRA's input Processors produce genes and proteins that carry both
UniProt and HGNC IDs in db_refs, whenever possible. FamPlex provides a name
space for protein families that are typically used in the literature.  More
information about FamPlex can be found here:
https://github.com/sorgerlab/famplex

+------------------------+------------------+--------------------------+
| Type                   | Database         | Example                  |
+========================+==================+==========================+
| Gene/Protein           | HGNC             | {'HGNC': '11998'}        |
+------------------------+------------------+--------------------------+
| Gene/Protein           | UniProt          | {'UP': 'P04637'}         |
+------------------------+------------------+--------------------------+
| Gene/Protein           | Entrez           | {'EGID': '5583'}         |
+------------------------+------------------+--------------------------+
| Gene/Protein family    | FamPlex          | {'FPLX': 'ERK'}          |
+------------------------+------------------+--------------------------+
| Gene/Protein family    | InterPro         | {'IP': 'IPR000308'}      |
+------------------------+------------------+--------------------------+
| Gene/Protein family    | Pfam             | {'PF': 'PF00071'}        |
+------------------------+------------------+--------------------------+
| Gene/Protein family    | NextProt family  | {'NXPFA': '03114'}       |
+------------------------+------------------+--------------------------+
| Chemical               | ChEBI            | {'CHEBI': 'CHEBI:63637'} |
+------------------------+------------------+--------------------------+
| Chemical               | PubChem          | {'PUBCHEM': '42611257'}  |
+------------------------+------------------+--------------------------+
| Chemical               | LINCS / HMS-LINCS| {'LINCS': '42611257'}    |
+------------------------+------------------+--------------------------+
| Metabolite             | HMDB             | {'HMDB': 'HMDB00122'}    |
+------------------------+------------------+--------------------------+
| Process, location, etc.| GO               | {'GO': 'GO:0006915'}     |
+------------------------+------------------+--------------------------+
| Process, disease, etc. | MeSH             | {'MESH': 'D008113'}      |
+------------------------+------------------+--------------------------+
| General terms          | NCIT             | {'NCIT': 'C28597'}       |
+------------------------+------------------+--------------------------+
| Raw text               | TEXT             | {'TEXT': 'Nf-kappaB'}    |
+------------------------+------------------+--------------------------+
"""

from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible

__all__ = [
    # Condition classes
    'BoundCondition', 'MutCondition', 'ModCondition', 'ActivityCondition',

    # Statement classes
    'Statement', 'Modification', 'AddModification', 'RemoveModification',
    'SelfModification', 'Phosphorylation', 'Autophosphorylation',
    'Transphosphorylation', 'Dephosphorylation', 'Hydroxylation',
    'Dehydroxylation', 'Sumoylation', 'Desumoylation', 'Acetylation',
    'Deacetylation', 'Glycosylation', 'Deglycosylation', 'Ribosylation',
    'Deribosylation', 'Ubiquitination', 'Deubiquitination', 'Farnesylation',
    'Defarnesylation', 'Geranylgeranylation', 'Degeranylgeranylation',
    'Palmitoylation', 'Depalmitoylation', 'Myristoylation', 'Demyristoylation',
    'Methylation', 'Demethylation', 'RegulateActivity', 'Inhibition',
    'Activation', 'GtpActivation', 'ActiveForm', 'HasActivity', 'Gef', 'Gap',
    'Complex', 'Translocation', 'RegulateAmount', 'DecreaseAmount',
    'IncreaseAmount', 'Conversion',

    # Error classes
    'InvalidLocationError',
    'InvalidResidueError', 'NotAStatementName',

    # Other classes
    'Agent',

    # Functions and values
    'get_valid_residue', 'get_valid_location', 'get_valid_location',
    'amino_acids', 'amino_acids_reverse', 'activity_types',
    'cellular_components', 'cellular_components_reverse', 'modtype_to_modclass',
    'modclass_to_modtype', 'modtype_conditions', 'modtype_to_inverse',
    'modclass_to_inverse', 'get_statement_by_name', 'stmt_sbo_map'
    ]

import abc
import logging
import networkx
import itertools
from collections import OrderedDict as _o
from ..util import *
from .. import Statement
from indra.statements.bio.resources import *


logger = logging.getLogger(__name__)


@python_2_unicode_compatible
class Modification(Statement):
    """Generic statement representing the modification of a protein.

    Parameters
    ----------
    enz : :py:class:`indra.statement.Agent`
        The enzyme involved in the modification.
    sub : :py:class:`indra.statement.Agent`
        The substrate of the modification.
    residue : str or None
        The amino acid residue being modified, or None if it is unknown or
        unspecified.
    position : str or None
        The position of the modified amino acid, or None if it is unknown or
        unspecified.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """
    _agent_order = ['enz', 'sub']

    def __init__(self, enz, sub, residue=None, position=None, evidence=None):
        super(Modification, self).__init__(evidence)
        self.enz = enz
        self.sub = sub
        self.residue = get_valid_residue(residue)
        if isinstance(position, int):
            self.position = str(position)
        else:
            self.position = position

    def matches_key(self):
        if self.enz is None:
            enz_key = None
        else:
            enz_key = self.enz.matches_key()
        key = (type(self), enz_key, self.sub.matches_key(),
               str(self.residue), str(self.position))
        return str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Modification has two agents in agent_list.")
        self.enz = agent_list[0]
        self.sub = agent_list[1]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if self.enz is None and other.enz is None:
            enz_refinement = True
        elif self.enz is None and other.enz is not None:
            enz_refinement = False
        elif self.enz is not None and other.enz is None:
            enz_refinement = True
        else:
            enz_refinement = self.enz.refinement_of(other.enz, hierarchies)
        sub_refinement = self.sub.refinement_of(other.sub, hierarchies)
        if not (enz_refinement and sub_refinement):
            return False
        # For this to be a refinement of the other, the modifications either
        # have to match or have this one be a subtype of the other; in
        # addition, the sites have to match, or this one has to have site
        # information and the other one not.
        residue_matches = (other.residue is None or
                           (self.residue == other.residue))
        position_matches = (other.position is None or
                            (self.position == other.position))
        return residue_matches and position_matches

    def equals(self, other):
        matches = super(Modification, self).equals(other)
        matches = (matches and (self.residue == other.residue)
                   and (self.position == other.position))
        return matches

    def contradicts(self, other, hierarchies):
        # If the modifications are not the opposite polarity of the
        # same subtype
        if not modclass_to_inverse[self.__class__] == other.__class__:
            return False
        # Skip all instances of not fully specified modifications
        agents = (self.enz, self.sub, other.enz, other.sub)
        if not all(a is not None for a in agents):
            return False
        # If the entities don't match, they can't be contradicting
        # Here we check pairs of agents at each "position" and
        # make sure they are the same or they are refinements of each other
        for self_agent, other_agent in zip(self.agent_list(),
                                           other.agent_list()):
            if not (self_agent.entity_matches(other_agent) or \
                    self_agent.refinement_of(other_agent, hierarchies) or \
                    other_agent.refinement_of(self_agent, hierarchies)):
                return False
        # At this point the entities definitely match so we need to
        # check the specific site that is being modified
        if self.residue == other.residue and self.position == other.position:
            return True
        else:
            return False

    def _get_mod_condition(self):
        """Return a ModCondition corresponding to this Modification."""
        mod_type = modclass_to_modtype[self.__class__]
        if isinstance(self, RemoveModification):
            mod_type = modtype_to_inverse[mod_type]
        mc = ModCondition(mod_type, self.residue, self.position, True)
        return mc

    def to_json(self, use_sbo=False):
        generic = super(Modification, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        if self.enz is not None:
            json_dict['enz'] = self.enz.to_json()
            if use_sbo:
                # enzymatic catalyst
                json_dict['enz']['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000460'
        if self.sub is not None:
            json_dict['sub'] = self.sub.to_json()
            if use_sbo:
                # substrate
                json_dict['sub']['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000015'
        if self.residue is not None:
            json_dict['residue'] = self.residue
        if self.position is not None:
            json_dict['position'] = self.position
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        enz = json_dict.get('enz')
        sub = json_dict.get('sub')
        residue = json_dict.get('residue')
        position = json_dict.get('position')
        evidence = json_dict.get('evidence', [])
        if enz:
            enz = Agent._from_json(enz)
        if sub:
            sub = Agent._from_json(sub)
        stmt = cls(enz, sub, residue, position)
        return stmt

    def __str__(self):
        res_str = (', %s' % self.residue) if self.residue is not None else ''
        pos_str = (', %s' % self.position) if self.position is not None else ''
        s = ("%s(%s, %s%s%s)" %
                 (type(self).__name__, self.enz, self.sub, res_str, pos_str))
        return s


class AddModification(Modification):
    pass


class RemoveModification(Modification):
    pass


@python_2_unicode_compatible
class SelfModification(Statement):
    """Generic statement representing the self-modification of a protein.

    Parameters
    ----------
    enz : :py:class:`indra.statement.Agent`
        The enzyme involved in the modification, which is also the substrate.
    residue : str or None
        The amino acid residue being modified, or None if it is unknown or
        unspecified.
    position : str or None
        The position of the modified amino acid, or None if it is unknown or
        unspecified.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """
    _agent_order = ['enz']

    def __init__(self, enz, residue=None, position=None, evidence=None):
        super(SelfModification, self).__init__(evidence)
        self.enz = enz
        self.residue = get_valid_residue(residue)
        if isinstance(position, int):
            self.position = str(position)
        else:
            self.position = position

    def __str__(self):
        res_str = (', %s' % self.residue) if self.residue is not None else ''
        pos_str = (', %s' % self.position) if self.position is not None else ''
        s = ("%s(%s%s%s)" %
             (type(self).__name__, self.enz, res_str, pos_str))
        return s

    def matches_key(self):
        key = (type(self), self.enz.matches_key(),
               str(self.residue), str(self.position))
        return str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("SelfModification has one agent.")
        self.enz = agent_list[0]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.enz.refinement_of(other.enz, hierarchies):
            return False
        # For this to be a refinement of the other, the modifications either
        # have to match or have this one be a subtype of the other; in
        # addition, the sites have to match, or this one has to have site
        # information and the other one not.
        residue_matches = (other.residue is None or
                           (self.residue == other.residue))
        position_matches = (other.position is None or
                            (self.position == other.position))
        return residue_matches and position_matches

    def equals(self, other):
        matches = super(SelfModification, self).equals(other)
        matches = (matches and self.residue == other.residue
                   and self.position == other.position)
        return matches

    def _get_mod_condition(self):
        """Return a ModCondition corresponding to this Modification."""
        mod_type = modclass_to_modtype[self.__class__]
        mc = ModCondition(mod_type, self.residue, self.position, True)
        return mc

    def to_json(self, use_sbo=False):
        generic = super(SelfModification, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        if self.enz is not None:
            json_dict['enz'] = self.enz.to_json()
            if use_sbo:
                # enzymatic catalyst
                json_dict['enz']['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000460'
        if self.residue is not None:
            json_dict['residue'] = self.residue
        if self.position is not None:
            json_dict['position'] = self.position
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        enz = json_dict.get('enz')
        residue = json_dict.get('residue')
        position = json_dict.get('position')
        if enz:
            enz = Agent._from_json(enz)
        stmt = cls(enz, residue, position)
        return stmt


class Phosphorylation(AddModification):
    """Phosphorylation modification.

    Examples
    --------
    MEK (MAP2K1) phosphorylates ERK (MAPK1) at threonine 185:

    >>> mek = Agent('MAP2K1')
    >>> erk = Agent('MAPK1')
    >>> phos = Phosphorylation(mek, erk, 'T', '185')
    """
    pass


class Autophosphorylation(SelfModification):
    """Intramolecular autophosphorylation, i.e., in *cis*.

    Examples
    --------
    p38 bound to TAB1 cis-autophosphorylates itself (see :pmid:`19155529`).

    >>> tab1 = Agent('TAB1')
    >>> p38_tab1 = Agent('P38', bound_conditions=[BoundCondition(tab1)])
    >>> autophos = Autophosphorylation(p38_tab1)
    """
    pass


class Transphosphorylation(SelfModification):
    """Autophosphorylation in *trans.*

    Transphosphorylation assumes that a kinase is already bound to a substrate
    (usually of the same molecular species), and phosphorylates it in an
    intra-molecular fashion. The enz property of the statement must have
    exactly one bound_conditions entry, and we assume that enz phosphorylates
    this molecule. The bound_neg property is ignored here.
    """
    pass


class Dephosphorylation(RemoveModification):
    """Dephosphorylation modification.

    Examples
    --------
    DUSP6 dephosphorylates ERK (MAPK1) at T185:

    >>> dusp6 = Agent('DUSP6')
    >>> erk = Agent('MAPK1')
    >>> dephos = Dephosphorylation(dusp6, erk, 'T', '185')
    """
    pass


class Hydroxylation(AddModification):
    """Hydroxylation modification."""
    pass


class Dehydroxylation(RemoveModification):
    """Dehydroxylation modification."""
    pass


class Sumoylation(AddModification):
    """Sumoylation modification."""
    pass


class Desumoylation(RemoveModification):
    """Desumoylation modification."""
    pass


class Acetylation(AddModification):
    """Acetylation modification."""
    pass


class Deacetylation(RemoveModification):
    """Deacetylation modification."""
    pass


class Glycosylation(AddModification):
    """Glycosylation modification."""
    pass


class Deglycosylation(RemoveModification):
    """Deglycosylation modification."""
    pass


class Ribosylation(AddModification):
    """Ribosylation modification."""
    pass


class Deribosylation(RemoveModification):
    """Deribosylation modification."""
    pass


class Ubiquitination(AddModification):
    """Ubiquitination modification."""
    pass


class Deubiquitination(RemoveModification):
    """Deubiquitination modification."""
    pass


class Farnesylation(AddModification):
    """Farnesylation modification."""
    pass


class Defarnesylation(RemoveModification):
    """Defarnesylation modification."""
    pass


class Geranylgeranylation(AddModification):
    """Geranylgeranylation modification."""
    pass


class Degeranylgeranylation(RemoveModification):
    """Degeranylgeranylation modification."""
    pass


class Palmitoylation(AddModification):
    """Palmitoylation modification."""
    pass


class Depalmitoylation(RemoveModification):
    """Depalmitoylation modification."""
    pass


class Myristoylation(AddModification):
    """Myristoylation modification."""
    pass


class Demyristoylation(RemoveModification):
    """Demyristoylation modification."""
    pass


class Methylation(AddModification):
    """Methylation modification."""
    pass


class Demethylation(RemoveModification):
    """Demethylation modification."""
    pass


@python_2_unicode_compatible
class RegulateActivity(Statement):
    """Regulation of activity.

    This class implements shared functionality of Activation and Inhibition
    statements and it should not be instantiated directly.
    """

    # The constructor here is an abstractmethod so that this class cannot
    # be directly instantiated.
    __metaclass__ = abc.ABCMeta

    _agent_order = ['subj', 'obj']

    @abc.abstractmethod
    def __init__(self):
        pass

    def __setstate__(self, state):
        if 'subj_activity' in state:
            logger.warning('Pickle file is out of date!')
        state.pop('subj_activity', None)
        self.__dict__.update(state)

    def matches_key(self):
        key = (type(self), self.subj.matches_key(),
               self.obj.matches_key(), str(self.obj_activity),
               str(self.is_activation))
        return str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("%s has two agents." % self.__class__.__name__)
        self.subj = agent_list[0]
        self.obj = agent_list[1]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        if self.is_activation != other.is_activation:
            return False
        if self.subj.refinement_of(other.subj, hierarchies) and \
           self.obj.refinement_of(other.obj, hierarchies):
            obj_act_match = (self.obj_activity == other.obj_activity) or \
                hierarchies['activity'].isa('INDRA_ACTIVITIES',
                                            self.obj_activity,
                                            'INDRA_ACTIVITIES',
                                            other.obj_activity)
            if obj_act_match:
                return True
            else:
                return False
        else:
            return False

    def contradicts(self, other, hierarchies):
        # If they aren't opposite classes, it's not a contradiction
        if {self.__class__, other.__class__} != {Activation, Inhibition}:
            return False

        # If they aren't opposite classes, it's not a contradiction
        if self.is_activation == other.is_activation:
            return False
        # Skip all instances of not fully specified statements
        agents = (self.subj, self.obj, other.subj, other.obj)
        if not all(a is not None for a in agents):
            return False
        # If the entities don't match, they can't be contradicting
        # Here we check pairs of agents at each "position" and
        # make sure they are the same or they are refinements of each other
        for self_agent, other_agent in zip(self.agent_list(),
                                           other.agent_list()):
            if not (self_agent.entity_matches(other_agent) or \
                    self_agent.refinement_of(other_agent, hierarchies) or \
                    other_agent.refinement_of(self_agent, hierarchies)):
                return False
        # Otherwise they are contradicting
        return True

    def to_json(self, use_sbo=False):
        generic = super(RegulateActivity, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        if self.subj is not None:
            json_dict['subj'] = self.subj.to_json()
            if use_sbo:
                if self.is_activation:
                    json_dict['subj']['sbo'] = \
                        'http://identifiers.org/sbo/SBO:0000459'  # stimulator
                else:
                    json_dict['subj']['sbo'] = \
                        'http://identifiers.org/sbo/SBO:0000020'  # inhibitor
        if self.obj is not None:
            json_dict['obj'] = self.obj.to_json()
            if use_sbo:
                if self.is_activation:
                    json_dict['obj']['sbo'] = \
                        'http://identifiers.org/sbo/SBO:0000643'  # stimulated
                else:
                    json_dict['obj']['sbo'] = \
                        'http://identifiers.org/sbo/SBO:0000642'  # inhibited
        if self.obj_activity is not None:
            json_dict['obj_activity'] = self.obj_activity
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        subj = json_dict.get('subj')
        obj = json_dict.get('obj')
        obj_activity = json_dict.get('obj_activity')
        if subj:
            subj = Agent._from_json(subj)
        if obj:
            obj = Agent._from_json(obj)
        stmt = cls(subj, obj, obj_activity)
        return stmt

    def __str__(self):
        obj_act_str = ', %s' % self.obj_activity if \
            self.obj_activity != 'activity' else ''
        s = ("%s(%s, %s%s)" %
             (type(self).__name__, self.subj,
              self.obj, obj_act_str))
        return s

    def __repr__(self):
        return self.__str__()

    def equals(self, other):
        matches = super(RegulateActivity, self).equals(other)
        matches = (matches and self.obj_activity == other.obj_activity
                   and self.is_activation == other.is_activation)
        return matches

    def _get_activity_condition(self):
        """Return ActivityCondition corresponding to this RegulateActivity."""
        return ActivityCondition(self.obj_activity, self.is_activation)


class Inhibition(RegulateActivity):
    """Indicates that a protein inhibits or deactivates another protein.

    This statement is intended to be used for physical interactions where the
    mechanism of inhibition is not explicitly specified, which is often the
    case for descriptions of mechanisms extracted from the literature.

    Parameters
    ----------
    subj : :py:class:`Agent`
        The agent responsible for the change in activity, i.e., the "upstream"
        node.
    obj : :py:class:`Agent`
        The agent whose activity is influenced by the subject, i.e., the
        "downstream" node.
    obj_activity : Optional[str]
        The activity of the obj Agent that is affected, e.g., its "kinase"
        activity.
    evidence : list of :py:class:`Evidence`
        Evidence objects in support of the modification.
    """
    def __init__(self, subj, obj, obj_activity='activity', evidence=None):
        super(RegulateActivity, self).__init__(evidence)
        self.subj = subj
        self.obj = obj
        if obj_activity not in activity_types:
            logger.warning('Invalid activity type: %s' % obj_activity)
        self.obj_activity = obj_activity
        self.is_activation = False


class Activation(RegulateActivity):
    """Indicates that a protein activates another protein.

    This statement is intended to be used for physical interactions where the
    mechanism of activation is not explicitly specified, which is often the
    case for descriptions of mechanisms extracted from the literature.

    Parameters
    ----------
    subj : :py:class:`Agent`
        The agent responsible for the change in activity, i.e., the "upstream"
        node.
    obj : :py:class:`Agent`
        The agent whose activity is influenced by the subject, i.e., the
        "downstream" node.
    obj_activity : Optional[str]
        The activity of the obj Agent that is affected, e.g., its "kinase"
        activity.
    evidence : list of :py:class:`Evidence`
        Evidence objects in support of the modification.

    Examples
    --------

    MEK (MAP2K1) activates the kinase activity of ERK (MAPK1):

    >>> mek = Agent('MAP2K1')
    >>> erk = Agent('MAPK1')
    >>> act = Activation(mek, erk, 'kinase')
    """
    def __init__(self, subj, obj, obj_activity='activity', evidence=None):
        super(RegulateActivity, self).__init__(evidence)
        self.subj = subj
        self.obj = obj
        if obj_activity not in activity_types:
            logger.warning('Invalid activity type: %s' % obj_activity)
        self.obj_activity = obj_activity
        self.is_activation = True


class GtpActivation(Activation):
    pass


@python_2_unicode_compatible
class ActiveForm(Statement):
    """Specifies conditions causing an Agent to be active or inactive.

    Types of conditions influencing a specific type of biochemical activity can
    include modifications, bound Agents, and mutations.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The Agent in a particular active or inactive state. The sets
        of ModConditions, BoundConditions, and MutConditions on the given
        Agent instance indicate the relevant conditions.
    activity : str
        The type of activity influenced by the given set of conditions,
        e.g., "kinase".
    is_active : bool
        Whether the conditions are activating (True) or inactivating (False).
    """
    _agent_order = ['agent']

    def __init__(self, agent, activity, is_active, evidence=None):
        super(ActiveForm, self).__init__(evidence)
        self.agent = agent
        if agent.activity is not None:
            logger.warning('Agent in ActiveForm should not have ' +
                           'ActivityConditions.')
            agent.activity = None
        if activity not in activity_types:
            logger.warning('Invalid activity type: %s' % activity)
        self.activity = activity
        self.is_active = is_active

    def matches_key(self):
        key = (type(self), self.agent.matches_key(),
               str(self.activity), str(self.is_active))
        return str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("ActiveForm has one agent.")
        self.agent = agent_list[0]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.agent.refinement_of(other.agent, hierarchies):
            return False

        # Make sure that the relationships and activities match
        if (self.is_active == other.is_active) and \
            (self.activity == other.activity
             or hierarchies['activity'].isa('INDRA_ACTIVITIES', self.activity,
                                            'INDRA_ACTIVITIES', other.activity)):
                return True
        else:
            return False

    def contradicts(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check that the polarity is constradicting up front
        # TODO: we could also check for cases where the polarities are
        # the same but some of the state conditions have opposite
        # polarities, for instance, if the presence/lack of the
        # same modification activates a given agent, that could be
        # considered a contradiction.
        if self.is_active == other.is_active:
            return False
        # Check that the activity types are the same
        # TODO: we could check for refinements here
        if self.activity != other.activity:
            return False
        # If the agents are exactly the same, this is a contradiction
        if self.agent.matches_key() == other.agent.matches_key():
            return True
        # Otherwise, if the two agents are related at the level of entities
        # and their state is exactly the same, then they contradict
        if self.agent.state_matches_key() == other.agent.state_matches_key():
            if self.agent.isa(other.agent, hierarchies) or \
                other.agent.isa(self.agent, hierarchies):
                return True
        return False

    def to_json(self, use_sbo=False):
        generic = super(ActiveForm, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        json_dict.update({'agent': self.agent.to_json(),
                          'activity': self.activity,
                          'is_active': self.is_active})
        if use_sbo:
            json_dict['agent']['sbo'] = \
                'http://identifiers.org/sbo/SBO:0000644'  # modified
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        agent = json_dict.get('agent')
        if agent:
            agent = Agent._from_json(agent)
        else:
            logger.error('ActiveForm statement missing agent')
            return None
        activity = json_dict.get('activity')
        is_active = json_dict.get('is_active')
        if activity is None:
            logger.warning('ActiveForm activity missing, defaulting ' +
                           'to `activity`')
            activity = 'activity'
        if is_active is None:
            logger.warning('ActiveForm is_active missing, defaulting ' +
                           'to True')
            is_active = True
        stmt = cls(agent, activity, is_active)
        return stmt

    def __str__(self):
        s = ("ActiveForm(%s, %s, %s)" %
             (self.agent, self.activity, self.is_active))
        return s

    def equals(self, other):
        matches = super(ActiveForm, self).equals(other)
        matches = (matches and self.activity == other.activity
                   and self.is_active == other.is_active)
        return matches


@python_2_unicode_compatible
class HasActivity(Statement):
    """States that an Agent has or doesn't have a given activity type.

    With this Statement, one cane express that a given protein is a kinase, or,
    for instance, that it is a transcription factor. It is also possible to
    construct negative statements with which one epxresses, for instance,
    that a given protein is not a kinase.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The Agent that that statement is about. Note that the detailed state
        of the Agent is not relevant for this type of statement.
    activity : str
        The type of activity, e.g., "kinase".
    has_activity : bool
        Whether the given Agent has the given activity (True) or
        not (False).
    """
    _agent_order = ['agent']

    def __init__(self, agent, activity, has_activity, evidence=None):
        super(HasActivity, self).__init__(evidence)
        if agent.activity is not None:
            logger.warning('Agent in HasActivity should not have ' +
                           'ActivityConditions.')
            agent.activity = None
        self.agent = agent
        if activity not in activity_types:
            logger.warning('Invalid activity type: %s' % activity)
        self.activity = activity
        self.has_activity = has_activity

    def matches_key(self):
        key = (type(self), self.agent.matches_key(),
               str(self.activity), str(self.has_activity))
        return str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 1:
            raise ValueError("HasActivity has one agent.")
        self.agent = agent_list[0]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if not self.agent.refinement_of(other.agent, hierarchies):
            return False

        # Make sure that the relationships and activities match
        if (self.has_activity == other.has_activity) and \
            (self.activity == other.activity
             or hierarchies['activity'].isa(self.activity, other.activity)):
                return True
        else:
            return False

    def __str__(self):
        s = ("HasActivity(%s, %s, %s)" %
             (self.agent, self.activity, self.has_activity))
        return s

    def equals(self, other):
        matches = super(HasActivity, self).equals(other)
        matches = (matches and self.activity == other.activity
                   and self.has_activity == other.has_activity)
        return matches


@python_2_unicode_compatible
class Gef(Statement):
    """Exchange of GTP for GDP on a small GTPase protein mediated by a GEF.

    Represents the generic process by which a guanosine exchange factor (GEF)
    catalyzes nucleotide exchange on a GTPase protein.

    Parameters
    ----------
    gef : :py:class:`Agent`
        The guanosine exchange factor.
    ras : :py:class:`Agent`
        The GTPase protein.

    Examples
    --------
    SOS1 catalyzes nucleotide exchange on KRAS:

    >>> sos = Agent('SOS1')
    >>> kras = Agent('KRAS')
    >>> gef = Gef(sos, kras)
    """
    _agent_order = ['gef', 'ras']

    def __init__(self, gef, ras, evidence=None):
        super(Gef, self).__init__(evidence)
        self.gef = gef
        self.ras = ras

    def matches_key(self):
        key = (type(self), self.gef.matches_key(),
               self.ras.matches_key())
        return str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Gef has two agents.")
        self.gef = agent_list[0]
        self.ras = agent_list[1]

    def __str__(self):
        s = "Gef(%s, %s)" % (self.gef.name, self.ras.name)
        return s

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check the GEF
        if self.gef.refinement_of(other.gef, hierarchies) and \
           self.ras.refinement_of(other.ras, hierarchies):
            return True
        else:
            return False

    def equals(self, other):
        matches = super(Gef, self).equals(other)
        return matches

    def to_json(self, use_sbo=False):
        generic = super(Gef, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        if self.gef is not None:
            json_dict['gef'] = self.gef.to_json()
            if use_sbo:
                json_dict['gef']['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000013'  # catalyst
        if self.ras is not None:
            json_dict['ras'] = self.ras.to_json()
            if use_sbo:
                json_dict['ras']['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000015'  # substrate
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        gef = json_dict.get('gef')
        ras = json_dict.get('ras')
        if gef:
            gef = Agent._from_json(gef)
        if ras:
            ras = Agent._from_json(ras)
        stmt = cls(gef, ras)
        return stmt


@python_2_unicode_compatible
class Gap(Statement):
    """Acceleration of a GTPase protein's GTP hydrolysis rate by a GAP.

    Represents the generic process by which a GTPase activating protein (GAP)
    catalyzes GTP hydrolysis by a particular small GTPase protein.

    Parameters
    ----------
    gap : :py:class:`Agent`
        The GTPase activating protein.
    ras : :py:class:`Agent`
        The GTPase protein.

    Examples
    --------
    RASA1 catalyzes GTP hydrolysis on KRAS:

    >>> rasa1 = Agent('RASA1')
    >>> kras = Agent('KRAS')
    >>> gap = Gap(rasa1, kras)
    """
    _agent_order = ['gap', 'ras']

    def __init__(self, gap, ras, evidence=None):
        super(Gap, self).__init__(evidence)
        self.gap = gap
        self.ras = ras

    def matches_key(self):
        key = (type(self), self.gap.matches_key(),
               self.ras.matches_key())
        return str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("Gap has two agents.")
        self.gap = agent_list[0]
        self.ras = agent_list[1]

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check the GAP
        if self.gap.refinement_of(other.gap, hierarchies) and \
           self.ras.refinement_of(other.ras, hierarchies):
            return True
        else:
            return False

    def __str__(self):
        s = "Gap(%s, %s)" % (self.gap.name, self.ras.name)
        return s

    def equals(self, other):
        matches = super(Gap, self).equals(other)
        return matches

    def to_json(self, use_sbo=False):
        generic = super(Gap, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        if self.gap is not None:
            json_dict['gap'] = self.gap.to_json()
            if use_sbo:
                json_dict['gap']['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000013'  # catalyst
        if self.ras is not None:
            json_dict['ras'] = self.ras.to_json()
            if use_sbo:
                json_dict['ras']['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000015'  # substrate
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        gap = json_dict.get('gap')
        ras = json_dict.get('ras')
        if gap:
            gap = Agent._from_json(gap)
        if ras:
            ras = Agent._from_json(ras)
        stmt = cls(gap, ras)
        return stmt


@python_2_unicode_compatible
class Complex(Statement):
    """A set of proteins observed to be in a complex.

    Parameters
    ----------
    members : list of :py:class:`Agent`
        The set of proteins in the complex.

    Examples
    --------
    BRAF is observed to be in a complex with RAF1:

    >>> braf = Agent('BRAF')
    >>> raf1 = Agent('RAF1')
    >>> cplx = Complex([braf, raf1])
    """
    _agent_order = ['members']

    def __init__(self, members, evidence=None):
        super(Complex, self).__init__(evidence)
        self.members = members

    def matches_key(self):
        members = sorted(self.members, key=lambda x: x.matches_key())
        key = (type(self), tuple(m.matches_key() for m in members))
        return str(key)

    def entities_match_key(self):
        key = tuple(a.entity_matches_key() if a is not None
                    else None for a in sorted(self.members,
                                              key=lambda x: x.matches_key()))
        return key

    def set_agent_list(self, agent_list):
        self.members = agent_list

    def __str__(self):
        s = '%s(%s)' % (type(self).__name__,
                        (', '.join([('%s' % m) for m in self.members])))
        return s

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Make sure the length of the members list is the same. Note that this
        # treats Complex([A, B, C]) as distinct from Complex([A, B]), rather
        # than as a refinement.
        if len(self.members) != len(other.members):
            return False

        def match_members(self_members, other_members):
            # First build a bipartite graph of refinement links
            G = networkx.Graph()
            for (self_idx, self_member), (other_idx, other_member) in \
                itertools.product(enumerate(self_members),
                                  enumerate(other_members)):
                if self_member.refinement_of(other_member, hierarchies):
                    G.add_edge('S%d' % self_idx, 'O%d' % other_idx)
            # Then find a maximal matching in the bipartite graph
            match = networkx.algorithms.max_weight_matching(G)
            # If every member has a pair, it is a valid refinement
            return len(match) == len(self_members)

        return match_members(self.members, other.members)

    def equals(self, other):
        matches = super(Complex, self).equals(other)
        return matches

    def to_json(self, use_sbo=False):
        generic = super(Complex, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        members = [m.to_json() for m in self.members]
        json_dict['members'] = members
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        members = json_dict.get('members')
        members = [Agent._from_json(m) for m in members]
        stmt = cls(members)
        return stmt


@python_2_unicode_compatible
class Translocation(Statement):
    """The translocation of a molecular agent from one location to another.

    Parameters
    ----------
    agent : :py:class:`Agent`
        The agent which translocates.
    from_location : Optional[str]
        The location from which the agent translocates. This must
        be a valid GO cellular component name (e.g. "cytoplasm")
        or ID (e.g. "GO:0005737").
    to_location : Optional[str]
        The location to which the agent translocates. This must
        be a valid GO cellular component name or ID.
    """
    _agent_order = ['agent']

    def __init__(self, agent, from_location=None, to_location=None,
                 evidence=None):
        super(Translocation, self).__init__(evidence)
        self.agent = agent
        self.from_location = get_valid_location(from_location)
        self.to_location = get_valid_location(to_location)

    def set_agent_list(self, agent_list):
        if(len(agent_list) != 1):
            raise ValueError("Translocation has 1 agent")
        self.agent = agent_list[0]

    def __str__(self):
        s = ("Translocation(%s, %s, %s)" %
             (self.agent, self.from_location, self.to_location))
        return s

    def refinement_of(self, other, hierarchies=None):
        # Make sure the statement types match
        if type(self) != type(other):
            return False
        # Check several conditions for refinement
        ch = hierarchies['cellular_component']
        ref1 = self.agent.refinement_of(other.agent, hierarchies)
        ref2 = (other.from_location is None or
                self.from_location == other.from_location or
                ch.partof('INDRA_LOCATIONS', self.from_location,
                          'INDRA_LOCATIONS', other.from_location))
        ref3 = (other.to_location is None or
                self.to_location == other.to_location or
                ch.partof('INDRA_LOCATIONS', self.to_location,
                          'INDRA_LOCATIONS', other.to_location))
        return (ref1 and ref2 and ref3)

    def equals(self, other):
        matches = super(Translocation, self).equals(other)
        matches = matches and (self.from_location == other.from_location)
        matches = matches and (self.to_location == other.to_location)
        return matches

    def matches_key(self):
        key = (type(self), self.agent.matches_key(), str(self.from_location),
               str(self.to_location))
        return str(key)

    def to_json(self, use_sbo=False):
        generic = super(Translocation, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        json_dict['agent'] = self.agent.to_json()
        if self.from_location is not None:
            json_dict['from_location'] = self.from_location
        if self.to_location is not None:
            json_dict['to_location'] = self.to_location
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        agent = json_dict.get('agent')
        if agent:
            agent = Agent._from_json(agent)
        else:
            logger.error('Translocation statement missing agent')
            return None
        from_location = json_dict.get('from_location')
        to_location = json_dict.get('to_location')
        stmt = cls(agent, from_location, to_location)
        return stmt


@python_2_unicode_compatible
class RegulateAmount(Statement):
    """Superclass handling operations on directed, two-element interactions."""
    _agent_order = ['subj', 'obj']

    def __init__(self, subj, obj, evidence=None):
        super(RegulateAmount, self).__init__(evidence)
        self.subj = subj
        if obj is None:
            raise ValueError('Object of %s cannot be None.' %
                             type(self).__name__)
        self.obj = obj

    def matches_key(self):
        if self.subj is None:
            subj_key = None
        else:
            subj_key = self.subj.matches_key()
        key = (type(self), subj_key, self.obj.matches_key())
        return str(key)

    def set_agent_list(self, agent_list):
        if len(agent_list) != 2:
            raise ValueError("%s has two agents in agent_list." %
                             type(self).__name__)
        self.subj = agent_list[0]
        self.obj = agent_list[1]

    def to_json(self, use_sbo=False):
        generic = super(RegulateAmount, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        if self.subj is not None:
            json_dict['subj'] = self.subj.to_json()
            if use_sbo:
                if isinstance(self, IncreaseAmount):
                    json_dict['subj']['sbo'] = \
                        'http://identifiers.org/sbo/SBO:0000459'  # stimulator
                else:
                    json_dict['subj']['sbo'] = \
                        'http://identifiers.org/sbo/SBO:0000020'  # inhibitor
        if self.obj is not None:
            json_dict['obj'] = self.obj.to_json()
            if use_sbo:
                if isinstance(self, IncreaseAmount):
                    json_dict['obj']['sbo'] = \
                        'http://identifiers.org/sbo/SBO:0000011'  # product
                else:
                    json_dict['obj']['sbo'] = \
                        'http://identifiers.org/sbo/SBO:0000010'  # reactant
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        subj = json_dict.get('subj')
        obj = json_dict.get('obj')
        if subj:
            subj = Agent._from_json(subj)
        if obj:
            obj = Agent._from_json(obj)
        stmt = cls(subj, obj)
        return stmt

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        if self.subj is None and other.subj is None:
            subj_refinement = True
        elif self.subj is None and other.subj is not None:
            subj_refinement = False
        elif self.subj is not None and other.subj is None:
            subj_refinement = True
        else:
            subj_refinement = self.subj.refinement_of(other.subj, hierarchies)
        obj_refinement = self.obj.refinement_of(other.obj, hierarchies)
        return (subj_refinement and obj_refinement)

    def equals(self, other):
        matches = super(RegulateAmount, self).equals(other)
        return matches

    def contradicts(self, other, hierarchies):
        # If they aren't opposite classes, it's not a contradiction
        if {self.__class__, other.__class__} != \
            {IncreaseAmount, DecreaseAmount}:
            return False
        # Skip all instances of not fully specified statements
        agents = (self.subj, self.obj, other.subj, other.obj)
        if not all(a is not None for a in agents):
            return False
        # If the entities don't match, they can't be contradicting
        # Here we check pairs of agents at each "position" and
        # make sure they are the same or they are refinements of each other
        for self_agent, other_agent in zip(self.agent_list(),
                                           other.agent_list()):
            if not (self_agent.entity_matches(other_agent) or \
                    self_agent.refinement_of(other_agent, hierarchies) or \
                    other_agent.refinement_of(self_agent, hierarchies)):
                return False
        # Otherwise they are contradicting
        return True

    def __str__(self):
        s = ("%s(%s, %s)" % (type(self).__name__, self.subj, self.obj))
        return s


class DecreaseAmount(RegulateAmount):
    """Degradation of a protein, possibly mediated by another protein.

    Note that this statement can also be used to represent inhibitors of
    synthesis (e.g., cycloheximide).

    Parameters
    ----------
    subj : :py:class:`indra.statement.Agent`
        The protein mediating the degradation.
    obj : :py:class:`indra.statement.Agent`
        The protein that is degraded.
    evidence : list of :py:class:`Evidence`
        Evidence objects in support of the degradation statement.
    """
    pass


class IncreaseAmount(RegulateAmount):
    """Synthesis of a protein, possibly mediated by another protein.

    Parameters
    ----------
    subj : :py:class:`indra.statement.Agent`
        The protein mediating the synthesis.
    obj : :py:class:`indra.statement.Agent`
        The protein that is synthesized.
    evidence : list of :py:class:`Evidence`
        Evidence objects in support of the synthesis statement.
    """
    pass


class Conversion(Statement):
    """Conversion of molecular species mediated by a controller protein.

    Parameters
    ----------
    subj : :py:class:`indra.statement.Agent`
        The protein mediating the conversion.
    obj_from : list of :py:class:`indra.statement.Agent`
        The list of molecular species being consumed by the conversion.
    obj_to : list of :py:class:`indra.statement.Agent`
        The list of molecular species being created by the conversion.
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        Evidence objects in support of the synthesis statement.
    """
    _agent_order = ['subj', 'obj_from', 'obj_to']

    def __init__(self, subj, obj_from=None, obj_to=None, evidence=None):
        super(Conversion, self).__init__(evidence=evidence)
        self.subj = subj
        self.obj_from = obj_from if obj_from is not None else []
        if isinstance(obj_from, Agent):
            self.obj_from = [obj_from]
        self.obj_to = obj_to if obj_to is not None else []
        if isinstance(obj_to, Agent):
            self.obj_to = [obj_to]

    def matches_key(self):
        keys = [type(self)]
        keys += [self.subj.matches_key() if self.subj else None]
        keys += [agent.matches_key() for agent in sorted_agents(self.obj_to)]
        keys += [agent.matches_key() for agent in sorted_agents(self.obj_from)]
        return str(keys)

    def set_agent_list(self, agent_list):
        num_obj_from = len(self.obj_from)
        num_obj_to = len(self.obj_to)
        if len(agent_list) != 1 + num_obj_from + num_obj_to:
            raise Exception('Conversion agent number must be preserved '
                            'when setting agent list.')
        self.subj = agent_list[0]
        self.obj_from = agent_list[1:num_obj_from+1]
        self.obj_to = agent_list[num_obj_from+1:]

    def to_json(self, use_sbo=False):
        generic = super(Conversion, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        if self.subj is not None:
            json_dict['subj'] = self.subj.to_json()
            if use_sbo:
                json_dict['subj']['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000013'  # catalyst
        json_dict['obj_from'] = [o.to_json() for o in self.obj_from]
        if use_sbo:
            for of in json_dict['obj_from']:
                of['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000010'  # reactant
        json_dict['obj_to'] = [o.to_json() for o in self.obj_to]
        if use_sbo:
            for ot in json_dict['obj_to']:
                ot['sbo'] = \
                    'http://identifiers.org/sbo/SBO:0000011'  # product
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        subj = json_dict.get('subj')
        obj_from = json_dict.get('obj_from')
        obj_to = json_dict.get('obj_to')
        if subj:
            subj = Agent._from_json(subj)
        if obj_from:
            obj_from = [Agent._from_json(o) for o in obj_from]
        if obj_to:
            obj_to = [Agent._from_json(o) for o in obj_to]
        stmt = cls(subj, obj_from, obj_to)
        return stmt

    def refinement_of(self, other, hierarchies):
        # Make sure the statement types match
        if type(self) != type(other):
            return False

        if self.subj is None and other.subj is None:
            subj_refinement = True
        elif self.subj is None and other.subj is not None:
            subj_refinement = False
        elif self.subj is not None and other.subj is None:
            subj_refinement = True
        else:
            subj_refinement = self.subj.refinement_of(other.subj, hierarchies)

        def refinement_agents(lst1, lst2):
            if len(lst1) != len(lst2):
                return False
            # Check that every agent in other is refined in self, but
            # only once!
            self_match_indices = set([])
            for other_agent in lst2:
                for self_agent_ix, self_agent in enumerate(lst1):
                    if self_agent_ix in self_match_indices:
                        continue
                    if self_agent.refinement_of(other_agent, hierarchies):
                        self_match_indices.add(self_agent_ix)
                        break
            if len(self_match_indices) != len(lst2):
                return False
            return True

        obj_from_refinement = refinement_agents(self.obj_from, other.obj_from)
        obj_to_refinement = refinement_agents(self.obj_to, other.obj_to)

        return (subj_refinement and obj_from_refinement and obj_to_refinement)

    def equals(self, other):
        matches = super(Conversion, self).equals(other)
        return matches

    def __str__(self):
        s = ("%s(%s, %s, %s)" % (type(self).__name__, self.subj, self.obj_from,
                                 self.obj_to))
        return s


# Mapping between modification type strings and subclasses of Modification
modtype_to_modclass = {str(cls.__name__.lower()): cls for cls in
                       AddModification.__subclasses__() +
                       RemoveModification.__subclasses__()}
# Add modification as a generic type
modtype_to_modclass['modification'] = Modification

modclass_to_modtype = {cls: str(cls.__name__.lower()) for cls in
                       AddModification.__subclasses__() +
                       RemoveModification.__subclasses__()}
# Add modification as a generic type
modclass_to_modtype[Modification] = 'modification'
modclass_to_modtype[Autophosphorylation] = 'phosphorylation'
modclass_to_modtype[Transphosphorylation] = 'phosphorylation'

# These are the modification types that are valid in ModConditions
modtype_conditions = [modclass_to_modtype[mt] for mt in
                      AddModification.__subclasses__()]
modtype_conditions.append('modification')


def _get_mod_inverse_maps():
    modtype_to_inverse = {}
    modclass_to_inverse = {}
    for cls in AddModification.__subclasses__():
        modtype = modclass_to_modtype[cls]
        modtype_inv = 'de' + modtype
        cls_inv = modtype_to_modclass[modtype_inv]
        modtype_to_inverse[modtype] = modtype_inv
        modtype_to_inverse[modtype_inv] = modtype
        modclass_to_inverse[cls] = cls_inv
        modclass_to_inverse[cls_inv] = cls
    return modtype_to_inverse, modclass_to_inverse


modtype_to_inverse, modclass_to_inverse = _get_mod_inverse_maps()


from indra.statements.bio.agent import *


stmt_sbo_map = {
    'acetylation': '0000215',
    'glycosylation': '0000217',
    'hydroxylation': '0000233',
    'methylation': '0000214',
    'myristoylation': '0000219',
    'palmitoylation': '0000218',
    'phosphorylation': '0000216',
    'farnesylation': '0000222',
    'geranylgeranylation': '0000223',
    'ubiquitination': '0000224',
    'dephosphorylation': '0000330',
    'addmodification': '0000210',  # addition of a chemical group
    'removemodification': '0000211',  # removal of a chemical group
    'modification': '0000182',  # conversion
    'conversion': '0000182',  # conversion
    'autophosphorylation': '0000216',  # phosphorylation
    'transphosphorylation': '0000216',  # phosphorylation
    'decreaseamount': '0000179',  # degradation
    'increaseamount': '0000183',  # transcription
    'complex': '0000526',  # protein complex formation
    'translocation': '0000185',  # transport reaction
    'regulateactivity': '0000182',  # conversion
    'activeform': '0000412',  # biological activity
    'rasgef': '0000172',  # catalysis
    'rasgap': '0000172',  # catalysis
    'statement': '0000231'  # occuring entity representation
    }
