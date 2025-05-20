"""
This is a minimal copy/modification of the tree related code from ete3.
Please see ete3 for more detail and the original code.
LB: This is here because I don't need the full dependency and also don't want to have to wait for
them to update to newer versions or fix problems...
"""
import os
import re
from collections import deque

DEFAULT_DIST = 1.0
DEFAULT_SUPPORT = 1.0
DEFAULT_NAME = ""

ITERABLE_TYPES = {list, set, tuple, frozenset}
_ILEGAL_NEWICK_CHARS = ":;(),[]\t\n\r="
_NHX_RE = "[&&NHX:[^]]*]"

_QUOTED_TEXT_RE = r"""((?=["'])(?:"[^"\\]*(?:\\[\s\S][^"\\]*)*"|'[^'\\]*(?:\\[\s\S][^'\\]*)*'))"""
_QUOTED_TEXT_PREFIX = 'ete3_quotref_'
_NAME_RE = "[^():,;]+?"
_FLOAT_RE = r"\s*[+-]?\d+\.?\d*(?:[eE][-+]?\d+)?\s*"
FLOAT_FORMATTER = "%0.6g"
NAME_FORMATTER = "%s"


class TreeError(Exception):
    """
    A problem occurred during a TreeNode operation
    """

    def __init__(self, value=''):
        self.value = value

    def __str__(self):
        return repr(self.value)


class TreeNode(object):
    def __init__(self, newick=None, format=0, dist=None, support=None,
                 name=None, quoted_node_names=False):
        self._children = []
        self._up = None
        self._dist = DEFAULT_DIST
        self._support = DEFAULT_SUPPORT
        self._img_style = None
        self.features = set([])
        # Add basic features
        self.features.update(["dist", "support", "name"])
        if dist is not None:
            self.dist = dist
        if support is not None:
            self.support = support

        self.name = name if name is not None else DEFAULT_NAME

        # Initialize tree
        if newick is not None:
            self._dist = 0.0
            read_newick(newick, root_node=self, format=format,
                        quoted_names=quoted_node_names)

    def __len__(self):
        """Node len returns number of children."""
        return len(self.get_leaves())

    def iter_leaves(self, is_leaf_fn=None):
        """
        Returns an iterator over the leaves under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        for n in self.traverse(strategy="preorder", is_leaf_fn=is_leaf_fn):
            if not is_leaf_fn:
                if n.is_leaf():
                    yield n
            else:
                if is_leaf_fn(n):
                    yield n

    def get_leaves(self, is_leaf_fn=None):
        """
        Returns the list of terminal nodes (leaves) under this node.

        :argument None is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.
        """
        return [n for n in self.iter_leaves(is_leaf_fn=is_leaf_fn)]

    def __iter__(self):
        """ Iterator over leaf nodes"""
        return self.iter_leaves()

    @property
    def children(self):
        return self._children

    @children.setter
    def children(self, value):
        if isinstance(value, list):
            for n in value:
                if not isinstance(n, type(self)):
                    raise TreeError("Incorrect child node type")
            self._children = value
        else:
            raise TreeError("Incorrect children type")

    @property
    def up(self):
        return self._up

    @up.setter
    def up(self, value):
        if isinstance(value, type(self)) or value is None:
            self._up = value
        else:
            raise TreeError("bad node_up type")

    @property
    def support(self):
        return self._support

    @support.setter
    def support(self, value):
        try:
            self._support = float(value)
        except (ValueError, TypeError):
            raise TreeError("node support must be a float number")

    # Topology management
    def add_child(self, child=None, name=None, dist=None, support=None):
        """
        Adds a new child to this node. If child node is not suplied
        as an argument, a new node instance will be created.

        :argument None child: the node instance to be added as a child.
        :argument None name: the name that will be given to the child.
        :argument None dist: the distance from the node to the child.
        :argument None support: the support value of child partition.

        :returns: The child node instance

        """
        if child is None:
            child = self.__class__()

        if name is not None:
            child.name = name
        if dist is not None:
            child.dist = dist
        if support is not None:
            child.support = support

        self.children.append(child)
        child.up = self
        return child

    def add_feature(self, pr_name, pr_value):
        """
        Add or update a node's feature.
        """
        setattr(self, pr_name, pr_value)
        self.features.add(pr_name)

    def add_features(self, **features):
        """
        Add or update several features.
        """
        for fname, fvalue in features.items():
            setattr(self, fname, fvalue)
            self.features.add(fname)

    def traverse(self, strategy="levelorder", is_leaf_fn=None):
        """
        Returns an iterator to traverse the tree structure under this
        node.

        :argument "levelorder" strategy: set the way in which tree
           will be traversed. Possible values are: "preorder" (first
           parent and then children) 'postorder' (first children and
           the parent) and "levelorder" (nodes are visited in order
           from root to leaves)

        :argument None is_leaf_fn: If supplied, ``is_leaf_fn``
           function will be used to interrogate nodes about if they
           are terminal or internal. ``is_leaf_fn`` function should
           receive a node instance as first argument and return True
           or False. Use this argument to traverse a tree by
           dynamically collapsing internal nodes matching
           ``is_leaf_fn``.
        """
        if strategy == "preorder":
            return self._iter_descendants_preorder(is_leaf_fn=is_leaf_fn)
        elif strategy == "levelorder":
            return self._iter_descendants_levelorder(is_leaf_fn=is_leaf_fn)
        elif strategy == "postorder":
            return self._iter_descendants_postorder(is_leaf_fn=is_leaf_fn)

    def iter_prepostorder(self, is_leaf_fn=None):
        """
        Iterate over all nodes in a tree yielding every node in both
        pre and post order. Each iteration returns a postorder flag
        (True if node is being visited in postorder) and a node
        instance.
        """
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                yield (False, node)
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
            else:
                # POSTORDER ACTIONS
                yield (True, node)

    def _iter_descendants_postorder(self, is_leaf_fn=None):
        to_visit = [self]
        if is_leaf_fn is not None:
            _leaf = is_leaf_fn
        else:
            _leaf = self.__class__.is_leaf

        while to_visit:
            node = to_visit.pop(-1)
            try:
                node = node[1]
            except TypeError:
                # PREORDER ACTIONS
                if not _leaf(node):
                    # ADD CHILDREN
                    to_visit.extend(reversed(node.children + [[1, node]]))
                else:
                    yield node
            else:
                # POSTORDER ACTIONS
                yield node

    def _iter_descendants_levelorder(self, is_leaf_fn=None):
        """
        Iterate over all desdecendant nodes.
        """
        tovisit = deque([self])
        while len(tovisit) > 0:
            node = tovisit.popleft()
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                tovisit.extend(node.children)

    def _iter_descendants_preorder(self, is_leaf_fn=None):
        """
        Iterator over all descendant nodes.
        """
        to_visit = deque()
        node = self
        while node is not None:
            yield node
            if not is_leaf_fn or not is_leaf_fn(node):
                to_visit.extendleft(reversed(node.children))
            try:
                node = to_visit.popleft()
            except:
                node = None

    def get_distance(self, target, target2=None, topology_only=False):
        """
        Returns the distance between two nodes. If only one target is
        specified, it returns the distance between the target and the
        current node.

        :argument target: a node within the same tree structure.

        :argument target2: a node within the same tree structure. If
          not specified, current node is used as target2.

        :argument False topology_only: If set to True, distance will
          refer to the number of nodes between target and target2.

        :returns: branch length distance between target and
          target2. If topology_only flag is True, returns the number
          of nodes between target and target2.

        """

        if target2 is None:
            target2 = self
            root = self.get_tree_root()
        else:
            # is target node under current node?
            root = self

        target, target2 = _translate_nodes(root, target, target2)
        ancestor = root.get_common_ancestor(target, target2)

        dist = 0.0
        for n in [target2, target]:
            current = n
            while current != ancestor:
                if topology_only:
                    dist += 1
                else:
                    dist += current.dist
                current = current.up
        if topology_only and target != target2:
            # counted ancestor once more than needed in while loop
            dist -= 1
        return dist

    def robinson_foulds(self, t2, attr_t1="name", attr_t2="name"):
        """
        Returns the Robinson-Foulds symmetric distance between current
        tree and a different tree instance.

        :param t2: reference tree

        :param name attr_t1: Compare trees using a custom node
                              attribute as a node name.

        :param name attr_t2: Compare trees using a custom node
                              attribute as a node name in target tree.

        :returns: (rf, rf_max, common_attrs, names, edges_t1, edges_t2)
        """
        ref_t = self
        target_t = t2
        if len(ref_t.children) > 2 or len(target_t.children) > 2:
            raise TreeError(
                "Unrooted tree found! Not supported in this package...")

        attrs_t1 = {getattr(n, attr_t1) for n in ref_t.iter_leaves() if hasattr(n, attr_t1)}
        attrs_t2 = {getattr(n, attr_t2) for n in target_t.iter_leaves() if hasattr(n, attr_t2)}
        if attrs_t1 != attrs_t2:
            raise TreeError("Trees have different taxa sets, not supported...")

        edges1 = extract_edge_set(ref_t, attr_t1)
        edges2 = extract_edge_set(target_t, attr_t2)

        # the two root edges are never counted here, as they are always
        # present in both trees because of the common attr filters
        rf = len((edges1 ^ edges2))

        # Otherwise we need to count the actual number of valid
        # partitions in each tree -2 is to avoid counting the root
        # partition of the two trees (only needed in rooted trees)
        max_parts = (sum(1 for p in edges1 if len(p) > 1) +
                     sum(1 for p in edges2 if len(p) > 1) - 2)

        min_comparison = (rf, max_parts, edges1, edges2)

        return min_comparison

    def get_cached_content(self, store_attr=None, container_type=set, leaves_only=True,
                           _store=None):
        """
        Returns a dictionary pointing to the preloaded content of each
        internal node under this tree. Such a dictionary is intended
        to work as a cache for operations that require many traversal
        operations.

        :param None store_attr: Specifies the node attribute that
            should be cached (i.e. name, distance, etc.). When none,
            the whole node instance is cached.

        :param _store: (internal use)
        """

        if _store is None:
            _store = {}

        def get_value(_n):
            if store_attr is None:
                _val = [_n]
            else:
                if not isinstance(store_attr, str):
                    _val = [tuple(getattr(_n, attr, None) for attr in store_attr)]

                else:
                    _val = [getattr(_n, store_attr, None)]

            return _val

        for ch in self.children:
            ch.get_cached_content(store_attr=store_attr,
                                  container_type=container_type,
                                  leaves_only=leaves_only,
                                  _store=_store)

        if self.children:
            if not leaves_only:
                val = container_type(get_value(self))
            else:
                val = container_type()
            for ch in self.children:
                if type(val) == list:
                    val.extend(_store[ch])
                if type(val) == set:
                    val.update(_store[ch])

                if not leaves_only:
                    if type(val) == list:
                        val.extend(get_value(ch))
                    if type(val) == set:
                        val.update(get_value(ch))

            _store[self] = val
        else:
            _store[self] = container_type(get_value(self))

        return _store

    def get_common_ancestor(self, *target_nodes, **kargs):
        """
        Returns the first common ancestor between this node and a given
        list of 'target_nodes'.

        **Examples:**

        ::

          t = tree.Tree("(((A:0.1, B:0.01):0.001, C:0.0001):1.0[&&NHX:name=common], (D:0.00001):0.000001):2.0[&&NHX:name=root];")
          A = t.get_descendants_by_name("A")[0]
          C = t.get_descendants_by_name("C")[0]
          common =  A.get_common_ancestor(C)
          print common.name

        """

        get_path = kargs.get("get_path", False)

        if len(target_nodes) == 1 and type(target_nodes[0]) \
                in set([set, tuple, list, frozenset]):
            target_nodes = target_nodes[0]

        # Convert node names into node instances
        target_nodes = _translate_nodes(self, *target_nodes)

        if type(target_nodes) != list:
            # If only one node is provided and is the same as the seed node,
            # return itself
            if target_nodes is self:
                if get_path:
                    return self, {}
                else:
                    return self
            else:
                # Otherwise find the common ancestor of current seed node and
                # the target_node provided
                target_nodes = [target_nodes, self]

        n2path = {}
        reference = []
        ref_node = None
        for n in target_nodes:
            current = n
            while current:
                n2path.setdefault(n, set()).add(current)
                if not ref_node:
                    reference.append(current)
                current = current.up
            if not ref_node:
                ref_node = n

        common = None
        for n in reference:
            broken = False
            for node, path in n2path.items():
                if node is not ref_node and n not in path:
                    broken = True
                    break

            if not broken:
                common = n
                break
        if not common:
            raise TreeError("Nodes are not connected!")

        if get_path:
            return common, n2path
        else:
            return common

    def get_tree_root(self):
        """
        Returns the absolute root node of current tree structure.
        """
        root = self
        while root.up is not None:
            root = root.up
        return root

    def is_leaf(self):
        """
        Return True if current node is a leaf.
        """
        return len(self.children) == 0

    def is_root(self):
        """
        Returns True if current node has no parent
        """
        if self.up is None:
            return True
        else:
            return False

    def write(self, features=None, outfile=None, format=0, is_leaf_fn=None,
              format_root_node=False, dist_formatter=None, support_formatter=None,
              name_formatter=None, quoted_node_names=False):
        """
        Returns the newick representation of current node. Several
        arguments control the way in which extra data is shown for
        every node:

        :argument features: a list of feature names to be exported
          using the Extended Newick Format (i.e. features=["name",
          "dist"]). Use an empty list to export all available features
          in each node (features=[])

        :argument outfile: writes the output to a given file

        :argument format: defines the newick standard used to encode the
          tree. See tutorial for details.

        :argument False format_root_node: If True, it allows features
          and branch information from root node to be exported as a
          part of the newick text string. For newick compatibility
          reasons, this is False by default.

        :argument is_leaf_fn: See :func:`TreeNode.traverse` for
          documentation.

        **Example:**

        ::

             t.write(features=["species","name"], format=1)

        """

        nw = write_newick(self, features=features, format=format,
                          is_leaf_fn=is_leaf_fn,
                          format_root_node=format_root_node,
                          dist_formatter=dist_formatter,
                          support_formatter=support_formatter,
                          name_formatter=name_formatter,
                          quoted_names=quoted_node_names)

        if outfile is not None:
            with open(outfile, "w") as OUT:
                OUT.write(nw)
        else:
            return nw


def extract_edge_set(t, attr_t):
    t_content = t.get_cached_content()
    edges = set()
    for content in t_content.values():
        edge = tuple(sorted(
            getattr(n, attr_t)
            for n in content
        ))
        if edge:
            edges.add(edge)
    return edges


def _translate_nodes(root, *nodes):
    name2node = dict([[n, None] for n in nodes if type(n) is str])
    if name2node:
        for n in root.traverse():
            if n.name in name2node:
                if name2node[n.name] is not None:
                    raise TreeError("Ambiguous node name: " + str(n.name))
                else:
                    name2node[n.name] = n

    if None in list(name2node.values()):
        notfound = [key for key, value in name2node.items() if value is None]
        raise ValueError("Node names not found: " + str(notfound))

    valid_nodes = []
    for n in nodes:
        if type(n) is not str:
            if type(n) is not root.__class__:
                raise TreeError("Invalid target node: " + str(n))
            else:
                valid_nodes.append(n)

    valid_nodes.extend(list(name2node.values()))
    if len(valid_nodes) == 1:
        return valid_nodes[0]
    else:
        return valid_nodes


class NewickError(Exception):
    """Exception class designed for NewickIO errors."""

    def __init__(self, value):
        if value is None:
            value = ''
        value += ("\nYou may want to check other newick loading flags"
                  " like 'format' or 'quoted_node_names'.")
        Exception.__init__(self, value)


def read_newick(newick, root_node=None, format=0, quoted_names=False):
    """ Reads a newick tree from either a string or a file, and returns
    an ETE tree structure.

    A previously existent node object can be passed as the root of the
    tree, which means that all its new children will belong to the same
    class as the root(This allows to work with custom TreeNode
    objects).

    You can also take advantage from this behaviour to concatenate
    several tree structures.
    """

    if root_node is None:
        root_node = TreeNode()

    if isinstance(newick, str):

        # try to determine whether the file exists.
        # For very large trees, if newick contains the content of the tree, rather than a file name,
        # this may fail, at least on Windows, because the os fails to stat the file content,
        # deeming it too long for testing with os.path.exists.  This raises a ValueError with
        # description "stat: path too long for Windows".  This is described in issue #258
        try:
            file_exists = os.path.exists(newick)
        except ValueError:  # failed to stat
            file_exists = False

        # if newick refers to a file, read it from file; otherwise,
        # regard it as a Newick content string.
        if file_exists:
            if newick.endswith('.gz'):
                import gzip
                with gzip.open(newick) as INPUT:
                    nw = INPUT.read()
            else:
                with open(newick) as INPUT:
                    nw = INPUT.read()
        else:
            nw = newick

        matcher = compile_matchers(formatcode=format)
        nw = nw.strip()
        if not nw.startswith('(') and nw.endswith(';'):
            # return _read_node_data(nw[:-1], root_node, "single", matcher, format)
            return _read_newick_from_string(nw, root_node, matcher, format, quoted_names)
        elif not nw.startswith('(') or not nw.endswith(';'):
            raise NewickError('Unexisting tree file or Malformed newick tree structure.')
        else:
            return _read_newick_from_string(nw, root_node, matcher, format, quoted_names)

    else:
        raise NewickError("'newick' argument must be either a filename or a newick string.")


NW_FORMAT = {
    0: [['name', str, True], ["dist", float, True], ['support', float, True],
        ["dist", float, True]],  # Flexible with support
    1: [['name', str, True], ["dist", float, True], ['name', str, True], ["dist", float, True]],
    # Flexible with internal node names
    2: [['name', str, False], ["dist", float, False], ['support', float, False],
        ["dist", float, False]],  # Strict with support values
    3: [['name', str, False], ["dist", float, False], ['name', str, False], ["dist", float, False]],
    # Strict with internal node names
    4: [['name', str, False], ["dist", float, False], [None, None, False], [None, None, False]],
    5: [['name', str, False], ["dist", float, False], [None, None, False], ["dist", float, False]],
    6: [['name', str, False], [None, None, False], [None, None, False], ["dist", float, False]],
    7: [['name', str, False], ["dist", float, False], ["name", str, False], [None, None, False]],
    8: [['name', str, False], [None, None, False], ["name", str, False], [None, None, False]],
    9: [['name', str, False], [None, None, False], [None, None, False], [None, None, False]],
    # Only topology with node names
    100: [[None, None, False], [None, None, False], [None, None, False], [None, None, False]]
    # Only Topology
}


def compile_matchers(formatcode):
    matchers = {}
    for node_type in ["leaf", "single", "internal"]:
        if node_type == "leaf" or node_type == "single":
            container1 = NW_FORMAT[formatcode][0][0]
            container2 = NW_FORMAT[formatcode][1][0]
            converterFn1 = NW_FORMAT[formatcode][0][1]
            converterFn2 = NW_FORMAT[formatcode][1][1]
            flexible1 = NW_FORMAT[formatcode][0][2]
            flexible2 = NW_FORMAT[formatcode][1][2]
        else:
            container1 = NW_FORMAT[formatcode][2][0]
            container2 = NW_FORMAT[formatcode][3][0]
            converterFn1 = NW_FORMAT[formatcode][2][1]
            converterFn2 = NW_FORMAT[formatcode][3][1]
            flexible1 = NW_FORMAT[formatcode][2][2]
            flexible2 = NW_FORMAT[formatcode][3][2]

        if converterFn1 == str:
            FIRST_MATCH = "(" + _NAME_RE + ")"
        elif converterFn1 == float:
            FIRST_MATCH = "(" + _FLOAT_RE + ")"
        elif converterFn1 is None:
            FIRST_MATCH = '()'

        if converterFn2 == str:
            SECOND_MATCH = "(:" + _NAME_RE + ")"
        elif converterFn2 == float:
            SECOND_MATCH = "(:" + _FLOAT_RE + ")"
        elif converterFn2 is None:
            SECOND_MATCH = '()'

        if flexible1 and node_type != 'leaf':
            FIRST_MATCH += "?"
        if flexible2:
            SECOND_MATCH += "?"

        matcher_str = r'^\s*%s\s*%s\s*(%s)?\s*$' % (FIRST_MATCH, SECOND_MATCH, _NHX_RE)
        compiled_matcher = re.compile(matcher_str)
        matchers[node_type] = [container1, container2, converterFn1, converterFn2, compiled_matcher]

    return matchers


def _read_newick_from_string(nw, root_node, matcher, formatcode, quoted_names):
    """ Reads a newick string in the New Hampshire format. """

    if quoted_names:
        # Quoted text is mapped to references
        quoted_map = {}
        unquoted_nw = ''
        counter = 0
        for token in re.split(_QUOTED_TEXT_RE, nw):
            counter += 1
            if counter % 2 == 1:  # normal newick tree structure data
                unquoted_nw += token
            else:  # quoted text, add to dictionary and replace with reference
                quoted_ref_id = _QUOTED_TEXT_PREFIX + str(int(counter / 2))
                unquoted_nw += quoted_ref_id
                quoted_map[quoted_ref_id] = token[1:-1]  # without the quotes
        nw = unquoted_nw

    if not nw.startswith('(') and nw.endswith(';'):
        _read_node_data(nw[:-1], root_node, "single", matcher, format)
        if quoted_names:
            if root_node.name.startswith(_QUOTED_TEXT_PREFIX):
                root_node.name = quoted_map[root_node.name]
        return root_node

    if nw.count('(') != nw.count(')'):
        raise NewickError('Parentheses do not match. Broken tree structure?')

    # white spaces and separators are removed
    nw = re.sub("[\n\r\t]+", "", nw)

    current_parent = None
    # Each chunk represents the content of a parent node, and it could contain
    # leaves and closing parentheses.
    # We may find:
    # leaf, ..., leaf,
    # leaf, ..., leaf))),
    # leaf)), leaf, leaf))
    # leaf))
    # ) only if formatcode == 100

    for chunk in nw.split("(")[1:]:
        # If no node has been created so far, this is the root, so use the node.
        current_parent = root_node if current_parent is None else current_parent.add_child()

        subchunks = [ch.strip() for ch in chunk.split(",")]
        # We should expect that the chunk finished with a comma (if next chunk
        # is an internal sister node) or a subchunk containing closing parenthesis
        # until the end of the tree.
        # [leaf, leaf, '']
        # [leaf, leaf, ')))', leaf, leaf, '']
        # [leaf, leaf, ')))', leaf, leaf, '']
        # [leaf, leaf, ')))', leaf), leaf, 'leaf);']
        if subchunks[-1] != '' and not subchunks[-1].endswith(';'):
            raise NewickError('Broken newick structure at: %s' % chunk)

        # lets process the subchunks. Every closing parenthesis will close a
        # node and go up one level.
        for i, leaf in enumerate(subchunks):
            if leaf.strip() == '' and i == len(subchunks) - 1:
                continue  # "blah blah ,( blah blah"
            closing_nodes = leaf.split(")")

            # first part after splitting by ) always contain leaf info
            _read_node_data(closing_nodes[0], current_parent, "leaf", matcher, formatcode)

            # next contain closing nodes and data about the internal nodes.
            if len(closing_nodes) > 1:
                for closing_internal in closing_nodes[1:]:
                    closing_internal = closing_internal.rstrip(";")
                    # read internal node data and go up one level
                    _read_node_data(closing_internal, current_parent, "internal", matcher,
                                    formatcode)
                    current_parent = current_parent.up

    # references in node names are replaced with quoted text before returning
    if quoted_names:
        for node in root_node.traverse():
            if node.name.startswith(_QUOTED_TEXT_PREFIX):
                node.name = quoted_map[node.name]

    return root_node


def _read_node_data(subnw, current_node, node_type, matcher, formatcode):
    """ Reads a leaf node from a subpart of the original newick
    tree """

    if node_type == "leaf" or node_type == "single":
        if node_type == "leaf":
            node = current_node.add_child()
        else:
            node = current_node
    else:
        node = current_node

    subnw = subnw.strip()

    if not subnw and node_type == 'leaf' and formatcode != 100:
        raise NewickError('Empty leaf node found')
    elif not subnw:
        return

    container1, container2, converterFn1, converterFn2, compiled_matcher = matcher[node_type]
    data = re.match(compiled_matcher, subnw)
    if data:
        data = data.groups()
        # This prevents ignoring errors even in flexible nodes:
        if subnw and data[0] is None and data[1] is None and data[2] is None:
            raise NewickError("Unexpected newick format '%s'" % subnw)

        if data[0] is not None and data[0] != '':
            node.add_feature(container1, converterFn1(data[0].strip()))

        if data[1] is not None and data[1] != '':
            node.add_feature(container2, converterFn2(data[1][1:].strip()))

        if data[2] is not None \
                and data[2].startswith("[&&NHX"):
            _parse_extra_features(node, data[2])
    else:
        raise NewickError("Unexpected newick format '%s' " % subnw[0:50])
    return


def _parse_extra_features(node, NHX_string):
    """ Reads node's extra data form its NHX string. NHX uses this
    format:  [&&NHX:prop1=value1:prop2=value2] """
    NHX_string = NHX_string.replace("[&&NHX:", "")
    NHX_string = NHX_string.replace("]", "")
    for field in NHX_string.split(":"):
        try:
            pname, pvalue = field.split("=")
        except ValueError as e:
            raise NewickError('Invalid NHX format %s' % field)
        node.add_feature(pname, pvalue)


def _read_node_data(subnw, current_node, node_type, matcher, formatcode):
    """ Reads a leaf node from a subpart of the original newick
    tree """

    if node_type == "leaf" or node_type == "single":
        if node_type == "leaf":
            node = current_node.add_child()
        else:
            node = current_node
    else:
        node = current_node

    subnw = subnw.strip()

    if not subnw and node_type == 'leaf' and formatcode != 100:
        raise NewickError('Empty leaf node found')
    elif not subnw:
        return

    container1, container2, converterFn1, converterFn2, compiled_matcher = matcher[node_type]
    data = re.match(compiled_matcher, subnw)
    if data:
        data = data.groups()
        # This prevents ignoring errors even in flexible nodes:
        if subnw and data[0] is None and data[1] is None and data[2] is None:
            raise NewickError("Unexpected newick format '%s'" % subnw)

        if data[0] is not None and data[0] != '':
            node.add_feature(container1, converterFn1(data[0].strip()))

        if data[1] is not None and data[1] != '':
            node.add_feature(container2, converterFn2(data[1][1:].strip()))

        if data[2] is not None \
                and data[2].startswith("[&&NHX"):
            _parse_extra_features(node, data[2])
    else:
        raise NewickError("Unexpected newick format '%s' " % subnw[0:50])
    return


def format_node(node, node_type, format, dist_formatter=None,
                support_formatter=None, name_formatter=None,
                quoted_names=False):
    if dist_formatter is None: dist_formatter = FLOAT_FORMATTER
    if support_formatter is None: support_formatter = FLOAT_FORMATTER
    if name_formatter is None: name_formatter = NAME_FORMATTER

    if node_type == "leaf":
        container1 = NW_FORMAT[format][0][0]  # name
        container2 = NW_FORMAT[format][1][0]  # dists
        converterFn1 = NW_FORMAT[format][0][1]
        converterFn2 = NW_FORMAT[format][1][1]
        flexible1 = NW_FORMAT[format][0][2]
    else:
        container1 = NW_FORMAT[format][2][0]  # support/name
        container2 = NW_FORMAT[format][3][0]  # dist
        converterFn1 = NW_FORMAT[format][2][1]
        converterFn2 = NW_FORMAT[format][3][1]
        flexible1 = NW_FORMAT[format][2][2]

    if converterFn1 == str:
        try:
            if not quoted_names:
                FIRST_PART = re.sub("[" + _ILEGAL_NEWICK_CHARS + "]", "_", \
                                    str(getattr(node, container1)))
            else:
                FIRST_PART = str(getattr(node, container1))
            if not FIRST_PART and container1 == 'name' and not flexible1:
                FIRST_PART = "NoName"

        except (AttributeError, TypeError):
            FIRST_PART = "?"

        FIRST_PART = name_formatter % FIRST_PART
        if quoted_names:
            # FIRST_PART = '"%s"' %FIRST_PART.decode('string_escape').replace('"', '\\"')
            FIRST_PART = '"%s"' % FIRST_PART

    elif converterFn1 is None:
        FIRST_PART = ""
    else:
        try:
            FIRST_PART = support_formatter % (converterFn2(getattr(node, container1)))
        except (ValueError, TypeError):
            FIRST_PART = "?"

    if converterFn2 == str:
        try:
            SECOND_PART = ":" + re.sub("[" + _ILEGAL_NEWICK_CHARS + "]", "_", \
                                       str(getattr(node, container2)))
        except (ValueError, TypeError):
            SECOND_PART = ":?"
    elif converterFn2 is None:
        SECOND_PART = ""
    else:
        try:
            # SECOND_PART = ":%0.6f" %(converterFn2(getattr(node, container2)))
            SECOND_PART = ":%s" % (dist_formatter % (converterFn2(getattr(node, container2))))
        except (ValueError, TypeError):
            SECOND_PART = ":?"

    return "%s%s" % (FIRST_PART, SECOND_PART)


def _get_features_string(self, features=None):
    """
    Generates the extended Newick NHX string with extra data about a node.
    """
    if features is None:
        features = []
    elif not features:
        features = sorted(self.features)

    parts = []
    for pr in features:
        if not hasattr(self, pr):
            continue

        raw = getattr(self, pr)

        if isinstance(raw, dict):
            raw = '|'.join(f"{k}-{v}" for k, v in raw.items())
        elif isinstance(raw, (list, tuple, set)):
            raw = '|'.join(map(str, raw))
        elif not isinstance(raw, str):
            raw = str(raw)

        # Replace illegal characters with "_"
        sanitized = re.sub(rf"[{_ILEGAL_NEWICK_CHARS}]", "_", raw)
        parts.append(f"{pr}={sanitized}")

    return f"[&{','.join(parts)}]" if parts else ""


def write_newick(rootnode, features=None, format=1, format_root_node=True,
                 is_leaf_fn=None, dist_formatter=None, support_formatter=None,
                 name_formatter=None, quoted_names=False):
    """
    Iteratively export a tree structure and returns its NHX
    representation.
    """
    newick = []
    leaf = is_leaf_fn if is_leaf_fn else lambda n: not bool(n.children)
    for postorder, node in rootnode.iter_prepostorder(is_leaf_fn=is_leaf_fn):
        if postorder:
            newick.append(")")
            if node.up is not None or format_root_node:
                newick.append(format_node(node, "internal", format,
                                          dist_formatter=dist_formatter,
                                          support_formatter=support_formatter,
                                          name_formatter=name_formatter,
                                          quoted_names=quoted_names))
                newick.append(_get_features_string(node, features))
        else:
            if node is not rootnode and node != node.up.children[0]:
                newick.append(",")

            if leaf(node):
                newick.append(format_node(node, "leaf", format,
                                          dist_formatter=dist_formatter,
                                          support_formatter=support_formatter,
                                          name_formatter=name_formatter,
                                          quoted_names=quoted_names))
                newick.append(_get_features_string(node, features))
            else:
                newick.append("(")

    newick.append(";")
    return ''.join(newick)


Tree = TreeNode
