from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible


__all__ = ['make_hash', 'sorted_agents', 'get_all_descendants',
           'get_type_hierarchy', 'NotAStatementName', 'get_statement_by_name',
           'make_statement_camel', 'get_unresolved_support_uuids']


from hashlib import md5
from . import Statement, Unresolved


def make_hash(s, n_bytes):
    """Make the hash from a matches key."""
    raw_h = int(md5(s.encode('utf-8')).hexdigest()[:n_bytes], 16)
    # Make it a signed int.
    return 16**n_bytes//2 - raw_h


def sorted_agents(agent_list):
    return sorted(agent_list, key=lambda ag: ag.matches_key())


def get_all_descendants(parent):
    """Get all the descendants of a parent class, recursively."""
    children = parent.__subclasses__()
    descendants = children[:]
    for child in children:
        descendants += get_all_descendants(child)
    return descendants


# In the future, when hierarchy is no longer determined by sub-classing, this
# function should be altered to account for the change.
def get_type_hierarchy(s):
    """Get the sequence of parents from `s` to Statement.

    Parameters
    ----------
    s : a class or instance of a child of Statement
        For example the statement `Phosphorylation(MEK(), ERK())` or just the
        class `Phosphorylation`.

    Returns
    -------
    parent_list : list[types]
        A list of the types leading up to Statement.

    Examples
    --------
        >> s = Phosphorylation(MAPK1(), Elk1())
        >> get_type_hierarchy(s)
        [Phosphorylation, AddModification, Modification, Statement]
        >> get_type_hierarchy(AddModification)
        [AddModification, Modification, Statement]
    """
    tp = type(s) if not isinstance(s, type) else s
    p_list = [tp]
    for p in tp.__bases__:
        if p is not Statement:
            p_list.extend(get_type_hierarchy(p))
        else:
            p_list.append(p)
    return p_list


class NotAStatementName(Exception):
    pass


def get_statement_by_name(stmt_name):
    """Get a statement class given the name of the statement class."""
    stmt_classes = get_all_descendants(Statement)
    for stmt_class in stmt_classes:
        if stmt_class.__name__.lower() == stmt_name.lower():
            return stmt_class
    raise NotAStatementName('\"%s\" is not recognized as a statement type!'
                            % stmt_name)


def make_statement_camel(stmt_name):
    """Makes a statement name match the case of the corresponding statement."""
    return get_statement_by_name(stmt_name).__name__


def get_unresolved_support_uuids(stmts):
    """Get uuids unresolved in support from stmts from stmts_from_json."""
    return {s.uuid for stmt in stmts for s in stmt.supports + stmt.supported_by
            if isinstance(s, Unresolved)}
