"""
Statements represent causal mechanisms between concepts of interest.

The evidence for a given Statement, which could include relevant citations,
database identifiers, and passages of text from the scientific literature, is
contained in one or more :py:class:`Evidence` objects associated with the
Statement.


JSON serialization of INDRA Statements
--------------------------------------

Statements can be serialized into JSON and deserialized from JSON to allow
their exchange in a platform-independent way. We also provide a JSON
schema (see http://json-schema.org to learn about schemas) in
https://raw.githubusercontent.com/sorgerlab/indra/master/indra/resources/statements_schema.json
which can be used to validate INDRA Statements JSONs.

Some validation tools include:

- jsonschema
a Python package to validate JSON content with respect to
a schema
- ajv-cli
Available at https://www.npmjs.com/package/ajv-cli
Install with "npm install -g ajv-cli" and then validate with:
    ajv -s statements_schema.json -d file_to_validate.json. This tool
provides more sophisticated and better interpretable output than
jsonschema.
- Web based tools
There are a variety of web-based tools for validation with JSON schemas,
including https://www.jsonschemavalidator.net
"""
import sys
import uuid
import networkx
from copy import deepcopy
from collections import OrderedDict as _o

try:  # Python 2
    basestring
except NameError:  # Python 3
    basestring = str


class Statement(object):
    """The parent class of all statements.

    Parameters
    ----------
    evidence : None or :py:class:`Evidence` or list of :py:class:`Evidence`
        If a list of Evidence objects is passed to the constructor, the
        value is set to this list. If a bare Evidence object is passed,
        it is enclosed in a list. If no evidence is passed (the default),
        the value is set to an empty list.
    supports : list of :py:class:`Statement`
        Statements that this Statement supports.
    supported_by : list of :py:class:`Statement`
        Statements supported by this statement.
    """
    _agent_order = NotImplemented

    def __init__(self, evidence=None, supports=None, supported_by=None):
        if evidence is None:
            self.evidence = []
        elif isinstance(evidence, Evidence):
            self.evidence = [evidence]
        elif isinstance(evidence, list):
            self.evidence = evidence
        else:
            raise ValueError('evidence must be an Evidence object, a list '
                             '(of Evidence objects), or None.')

        # Initialize supports/supported_by fields, which should be lists
        self.supports = supports if supports else []
        self.supported_by = supported_by if supported_by else []
        self.belief = 1
        self.uuid = '%s' % uuid.uuid4()
        self._full_hash = None
        self._shallow_hash = None
        return

    def matches_key(self):
        raise NotImplementedError("Method must be implemented in child class.")

    def matches(self, other):
        return self.matches_key() == other.matches_key()

    def get_hash(self, shallow=True, refresh=False):
        """Get a hash for this Statement.

        There are two types of hash, "shallow" and "full". A shallow hash is
        as unique as the information carried by the statement, i.e. it is a hash
        of the `matches_key`. This means that differences in source, evidence,
        and so on are not included. As such, it is a shorter hash (14 nibbles).
        The odds of a collision among all the statements we expect to encounter
        (well under 10^8) is ~10^-9 (1 in a billion). Checks for collisions can
        be done by using the matches keys.

        A full hash includes, in addition to the matches key, information from
        the evidence of the statement. These hashes will be equal if the two
        Statements came from the same sentences, extracted by the same reader,
        from the same source. These hashes are correspondingly longer (16
        nibbles). The odds of a collision for an expected less than 10^10
        extractions is ~10^-9 (1 in a billion).

        Note that a hash of the Python object will also include the `uuid`, so
        it will always be unique for every object.

        Parameters
        ----------
        shallow : bool
            Choose between the shallow and full hashes described above. Default
            is true (e.g. a shallow hash).
        refresh : bool
            Used to get a new copy of the hash. Default is false, so the hash,
            if it has been already created, will be read from the attribute.
            This is primarily used for speed testing.

        Returns
        -------
        hash : int
            A long integer hash.
        """
        if shallow:
            if not hasattr(self, '_shallow_hash') or self._shallow_hash is None \
                    or refresh:
                self._shallow_hash = make_hash(self.matches_key(), 14)
            ret = self._shallow_hash
        else:
            if not hasattr(self, '_full_hash') or self._full_hash is None \
                    or refresh:
                ev_mk_list = sorted([ev.matches_key() for ev in self.evidence])
                self._full_hash = \
                    make_hash(self.matches_key() + str(ev_mk_list), 16)
            ret = self._full_hash
        return ret

    def _tag_evidence(self):
        """Set all the Evidence stmt_tag to my deep matches-key hash."""
        h = self.get_hash(shallow=False)
        for ev in self.evidence:
            ev.stmt_tag = h
        return

    def agent_list_with_bound_condition_agents(self):
        # Returns the list of agents both directly participating in the
        # statement and referenced through bound conditions.
        l = self.agent_list()
        for a in self.agent_list():
            if a is not None:
                bc_agents = [bc.agent for bc in a.bound_conditions]
                l.extend(bc_agents)
        return l

    def agent_list(self, deep_sorted=False):
        """Get the canonicallized agent list."""
        ag_list = []
        for ag_name in self._agent_order:
            ag_attr = getattr(self, ag_name)
            if isinstance(ag_attr, Concept) or ag_attr is None:
                ag_list.append(ag_attr)
            elif isinstance(ag_attr, list):
                if not all([isinstance(ag, Concept) for ag in ag_attr]):
                    raise TypeError("Expected all elements of list to be Agent "
                                    "and/or Concept, but got: %s"
                                    % {type(ag) for ag in ag_attr})
                if deep_sorted:
                    ag_attr = sorted_agents(ag_attr)
                ag_list.extend(ag_attr)
            else:
                raise TypeError("Expected type Agent, Concept, or list, got "
                                "type %s." % type(ag_attr))
        return ag_list

    def entities_match(self, other):
        self_key = self.entities_match_key()
        other_key = other.entities_match_key()
        if len(self_key) != len(other_key):
            return False
        for self_agent, other_agent in zip(self_key, other_key):
            if self_agent is None or other_agent is None:
                continue
            if self_agent != other_agent:
                return False
        return True

    def entities_match_key(self):
        key = tuple(a.entity_matches_key() if a is not None
                    else None for a in self.agent_list())
        return key

    def print_supports(self):
        print('%s supported_by:' % str(self))
        if self.supported_by:
            print('-->')
            for s in self.supported_by:
                s.print_supports()

    def __repr__(self):
        if sys.version_info[0] >= 3:
            return str(self)
        else:
            return str(self).encode('utf-8')

    def equals(self, other):
        if type(self) != type(other):
            return False
        if len(self.agent_list()) == len(other.agent_list()):
            for s, o in zip(self.agent_list(), other.agent_list()):
                if (s is None and o is not None) or \
                        (s is not None and o is None):
                    return False
                if s is not None and o is not None and not s.equals(o):
                    return False
        else:
            return False
        if len(self.evidence) == len(other.evidence):
            for s, o in zip(self.evidence, other.evidence):
                if not s.equals(o):
                    return False
        else:
            return False
        return True

    def contradicts(self, other, hierarchies):
        # Placeholder for implementation in subclasses
        return False

    def to_json(self, use_sbo=False):
        """Return serialized Statement as a JSON dict.

        Parameters
        ----------
        use_sbo : Optional[bool]
            If True, SBO annotations are added to each applicable element of
            the JSON. Default: False

        Returns
        -------
        json_dict : dict
            The JSON-serialized INDRA Statement.
        """
        stmt_type = type(self).__name__
        # Original comment: For backwards compatibility, could be removed later
        all_stmts = [self] + self.supports + self.supported_by
        for st in all_stmts:
            if not hasattr(st, 'uuid'):
                st.uuid = '%s' % uuid.uuid4()
        ##################
        json_dict = _o({'type': stmt_type})
        json_dict['belief'] = self.belief
        if self.evidence:
            evidence = [ev.to_json() for ev in self.evidence]
            json_dict['evidence'] = evidence
        json_dict['id'] = '%s' % self.uuid
        if self.supports:
            json_dict['supports'] = \
                ['%s' % st.uuid for st in self.supports]
        if self.supported_by:
            json_dict['supported_by'] = \
                ['%s' % st.uuid for st in self.supported_by]

        def get_sbo_term(cls):
            sbo_term = stmt_sbo_map.get(cls.__name__.lower())
            while not sbo_term:
                cls = cls.__bases__[0]
                sbo_term = stmt_sbo_map.get(cls.__name__.lower())
            return sbo_term

        if use_sbo:
            sbo_term = get_sbo_term(self.__class__)
            json_dict['sbo'] = \
                'http://identifiers.org/sbo/SBO:%s' % sbo_term
        return json_dict

    @classmethod
    def _from_json(cls, json_dict):
        stmt_type = json_dict.get('type')
        stmt_cls = getattr(sys.modules[__name__], stmt_type)
        stmt = stmt_cls._from_json(json_dict)
        evidence = json_dict.get('evidence', [])
        stmt.evidence = [Evidence._from_json(ev) for ev in evidence]
        stmt.supports = json_dict.get('supports', [])[:]
        stmt.supported_by = json_dict.get('supported_by', [])[:]
        stmt.belief = json_dict.get('belief', 1.0)
        stmt_id = json_dict.get('id')
        if not stmt_id:
            stmt_id = '%s' % uuid.uuid4()
        stmt.uuid = stmt_id
        return stmt

    def to_graph(self):
        """Return Statement as a networkx graph."""
        def json_node(graph, element, prefix):
            if not element:
                return None
            node_id = '|'.join(prefix)
            if isinstance(element, list):
                graph.add_node(node_id, label='')
                # Enumerate children and add nodes and connect to anchor node
                for i, sub_element in enumerate(element):
                    sub_id = json_node(graph, sub_element, prefix + ['%s' % i])
                    if sub_id:
                        graph.add_edge(node_id, sub_id, label='')
            elif isinstance(element, dict):
                graph.add_node(node_id, label='')
                # Add node recursively for each element
                # Connect to this node with edge label according to key
                for k, v in element.items():
                    if k == 'id':
                        continue
                    elif k == 'name':
                        graph.node[node_id]['label'] = v
                        continue
                    elif k == 'type':
                        graph.node[node_id]['label'] = v
                        continue

                    sub_id = json_node(graph, v, prefix + ['%s' % k])
                    if sub_id:
                        graph.add_edge(node_id, sub_id, label=('%s' % k))
            else:
                if isinstance(element, basestring) and \
                        element.startswith('http'):
                    element = element.split('/')[-1]
                graph.add_node(node_id, label=('%s' % str(element)))
            return node_id
        jd = self.to_json()
        graph = networkx.DiGraph()
        json_node(graph, jd, ['%s' % self.uuid])
        return graph

    def make_generic_copy(self, deeply=False):
        """Make a new matching Statement with no provenance.

        All agents and other attributes besides evidence, belief, supports, and
        supported_by will be copied over, and a new uuid will be assigned.
        Thus, the new Statement will satisfy `new_stmt.matches(old_stmt)`.

        If `deeply` is set to True, all the attributes will be deep-copied,
        which is comparatively slow. Otherwise, attributes of this statement
        may be altered by changes to the new matching statement.
        """
        if deeply:
            kwargs = deepcopy(self.__dict__)
        else:
            kwargs = self.__dict__.copy()
        for attr in ['evidence', 'belief', 'uuid', 'supports', 'supported_by',
                     'is_activation']:
            kwargs.pop(attr, None)
        for attr in ['_full_hash', '_shallow_hash']:
            my_hash = kwargs.pop(attr, None)
            my_shallow_hash = kwargs.pop(attr, None)
        for attr in self._agent_order:
            attr_value = kwargs.get(attr)
            if isinstance(attr_value, list):
                kwargs[attr] = sorted_agents(attr_value)
        new_instance = self.__class__(**kwargs)
        new_instance._full_hash = my_hash
        new_instance._shallow_hash = my_shallow_hash
        return new_instance


class Unresolved(Statement):
    """A special statement type used in support when a uuid can't be resolved.

    When using the `stmts_from_json` method, it is sometimes not possible to
    resolve the uuid found in `support` and `supported_by` in the json
    representation of an indra statement. When this happens, this class is used
    as a place-holder, carrying only the uuid of the statement.
    """
    def __init__(self, uuid_str=None, shallow_hash=None, full_hash=None):
        super(Unresolved, self).__init__()
        self.uuid = uuid_str
        self._shallow_hash = shallow_hash
        self._full_hash = full_hash
        assert self.uuid or self._shallow_hash or self._full_hash, \
            "Some identifying information must be given."

    def __str__(self):
        if self.uuid:
            return "%s(uuid=%s)" % (type(self).__name__, self.uuid)
        elif self._shallow_hash:
            return "%s(shallow_hash=%s)" % (type(self).__name__,
                                            self._shallow_hash)
        else:
            return "%s(full_hash=%s)" % (type(self).__name__,
                                         self._full_hash)


from .io import *
from .util import *
from .context import *
from .evidence import *
from .bio.statements import *
from .general.statements import *


