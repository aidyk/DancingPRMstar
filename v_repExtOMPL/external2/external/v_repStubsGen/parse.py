from lxml import etree
import model

def parse(xml_file):
    tree = etree.parse(xml_file)
    root = tree.getroot()
    return model.Plugin(root)

