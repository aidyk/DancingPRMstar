class Enum(object):
    def __init__(self, plugin, node):
        if node.tag != 'enum':
            raise ValueError('expected <enum>, got <%s>' % node.tag)
        self.plugin = plugin
        self.name = node.attrib['name']
        self.item_prefix = node.attrib.get('item-prefix', '')
        self.base = node.attrib.get('base', None)
        self.items = [n.attrib['name'] for n in node.findall('item')]

