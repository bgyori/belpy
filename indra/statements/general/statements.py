from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str


__all__ = ['Influence', 'Association', 'Concept']

import sys
from collections import OrderedDict as _o
from indra.statements.general.concept import Concept
from indra.statements.bio.statements import IncreaseAmount, Complex


class Influence(IncreaseAmount):
    """An influence on the quantity of a concept of interest.

    Parameters
    ----------
    subj : :py:class:`indra.statement.Concept`
        The concept which acts as the influencer.
    obj : :py:class:`indra.statement.Concept`
        The concept which acts as the influencee
    subj_delta : Optional[dict]
        A dictionary specifying the polarity and magnitude of
        change in the subject.
    obj_delta : Optional[dict]
        A dictionary specifying the polarity and magnitude of
        change in the object.
    evidence : list of :py:class:`Evidence`
        Evidence objects in support of the statement.
    """
    def __init__(self, subj, obj, subj_delta=None, obj_delta=None,
                 evidence=None):
        super(Influence, self).__init__(subj, obj, evidence)
        if subj_delta is None:
            subj_delta = {'polarity': None, 'adjectives': []}
        if obj_delta is None:
            obj_delta = {'polarity': None, 'adjectives': []}
        self.subj_delta = subj_delta
        self.obj_delta = obj_delta

    def refinement_of(self, other, hierarchies):
        def delta_refinement(dself, dother):
            # Polarities are either equal
            if dself['polarity'] == dother['polarity']:
                pol_refinement = True
            # Or this one has a polarity and the other doesn't
            elif dself['polarity'] is not None and dother['polarity'] is None:
                pol_refinement = True
            else:
                pol_refinement = False

            # If other's adjectives are a subset of this
            if set(dother['adjectives']).issubset(set(dself['adjectives'])):
                adj_refinement = True
            else:
                adj_refinement = False
            return pol_refinement and adj_refinement

        # Make sure the statement types match
        if type(self) != type(other):
            return False

        # Check agent arguments
        subj_refinement = self.subj.refinement_of(other.subj, hierarchies)
        obj_refinement = self.obj.refinement_of(other.obj, hierarchies)
        subjd_refinement = delta_refinement(self.subj_delta, other.subj_delta)
        objd_refinement = delta_refinement(self.obj_delta, other.obj_delta)
        return (subj_refinement and obj_refinement and
                subjd_refinement and objd_refinement)

    def equals(self, other):
        def delta_equals(dself, dother):
            if (dself['polarity'] == dother['polarity']) and \
                    (set(dself['adjectives']) == set(dother['adjectives'])):
                return True
            else:
                return False
        matches = super(Influence, self).equals(other) and \
                  delta_equals(self.subj_delta, other.subj_delta) and \
                  delta_equals(self.obj_delta, other.obj_delta)
        return matches

    def matches_key(self):
        key = (type(self), self.subj.matches_key(),
               self.obj.matches_key(),
               self.subj_delta['polarity'],
               sorted(list(set(self.subj_delta['adjectives']))),
               self.obj_delta['polarity'],
               sorted(list(set(self.obj_delta['adjectives']))))
        return str(key)

    def contradicts(self, other, hierarchies):
        # First case is if they are "consistent" and related
        if self.entities_match(other) or \
                self.refinement_of(other, hierarchies) or \
                other.refinement_of(self, hierarchies):
            sp = self.overall_polarity()
            op = other.overall_polarity()
            if sp and op and sp * op == -1:
                return True
        # Second case is if they are "opposites" and related
        if (self.subj.entity_matches(other.subj) and \
            self.obj.is_opposite(other.obj, hierarchies)) or \
                (self.obj.entity_matches(other.obj) and \
                 self.subj.is_opposite(other.subj, hierarchies)):
            sp = self.overall_polarity()
            op = other.overall_polarity()
            if sp and op and sp * op == 1:
                return True
        return False

    def overall_polarity(self):
        # Set p1 and p2 to None / 1 / -1 depending on polarity
        p1 = self.subj_delta['polarity']
        p2 = self.obj_delta['polarity']
        if p1 is None and p2 is None:
            pol = None
        elif p2 is None:
            pol = p1
        elif p1 is None:
            pol = p2
        else:
            pol = p1 * p2
        return pol

    def to_json(self, use_sbo=False):
        generic = super(Influence, self).to_json(use_sbo)
        json_dict = _o({'type': generic['type']})
        json_dict['subj'] = generic['subj']
        json_dict['subj_delta'] = self.subj_delta
        json_dict['obj'] = generic['obj']
        json_dict['obj_delta'] = self.obj_delta
        json_dict.update(generic)
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        subj = json_dict.get('subj')
        obj = json_dict.get('obj')
        subj_delta = json_dict.get('subj_delta')
        obj_delta = json_dict.get('obj_delta')
        if subj:
            subj = Concept._from_json(subj)
        if obj:
            obj = Concept._from_json(obj)
        stmt = cls(subj, obj, subj_delta, obj_delta)
        return stmt

    def __repr__(self):
        if sys.version_info[0] >= 3:
            return self.__str__()
        else:
            return self.__str__().encode('utf-8')

    def __str__(self):
        def _influence_concept_str(concept, delta):
            if delta is not None:
                pol = delta.get('polarity')
                if pol == 1:
                    pol_str = 'positive'
                elif pol == -1:
                    pol_str = 'negative'
                else:
                    pol_str = ''
                concept_str = '%s(%s)' % (concept.name, pol_str)
            else:
                concept_str = concept.name
            return concept_str
        s = ("%s(%s, %s)" % (type(self).__name__,
                             _influence_concept_str(self.subj,
                                                    self.subj_delta),
                             _influence_concept_str(self.obj,
                                                    self.obj_delta)))
        return s


class Association(Complex):
    pass

