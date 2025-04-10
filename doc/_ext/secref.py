"""Custom Sphinx role. A mix between ref and numref.

Created 2023-08-01

@author: Michele Riva
"""

from sphinx.builders.latex import LaTeXBuilder


def secref_role(name, rawtext, text, lineno, inliner,
                options=None, content=None):
    """Return a numref in LaTeX, a ref otherwise.

    Parameters
    ----------
    name : str
        The role name used in the document.
    rawtext : str
        The entire role token, including the role name.
    text : str
        The parsed text in the token.
    lineno : int
        The line number of the role token.
    inliner : docutils.parsers.rst.states.Inliner
        The Sphinx inliner that calls this role resolver.
    options : dict, or None, optional
        Options for constructing the node. Default is None.
    content : list or None, optional
        The directive content for customization. Default is None.

    Returns
    -------
    nodes : list
        The new nodes created.
    error_messages : list
        Errors encountered while creating the nodes.
    """
    if options is None:
        options = {}
    if content is None:
        content = []

    env = inliner.document.settings.env
    std_domain = env.domains['std']

    # Pick a concrete reference role depending on the builder that is
    # currently in use. Also, stick to a simple ref in case there is
    # an explicit title given.
    builder = env.app.builder
    use_numref = (isinstance(builder, LaTeXBuilder)
                  and "<" not in text)
    new_name = "numref" if use_numref else "ref"

    ref_role = std_domain.roles[new_name]
    roles, messages = ref_role(f"std:{new_name}", rawtext, text,
                               lineno, inliner, options, content)

    # TODO: here we could also be able to preprocess any math that may
    # be present in the section header to which this node refers to.
    # roles[0] is a sphinx.addnodes.pending_xref, which is essentially
    # a docutils.nodes.Element. We may be able to replace some
    # attribute with the correctly processed return of a :math: role

    return roles, messages


def setup(app):
    """Add new role to Sphinx app."""
    app.add_role('secref', secref_role)
    return {
        'version': '1.0',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
        }
