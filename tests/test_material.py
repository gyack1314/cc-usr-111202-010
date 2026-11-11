"""核对来源读取、资料边界及原始摘要。"""
import json
import tempfile
import unittest
from pathlib import Path
from research_data.io import finite,integer,rational,load,unique
from research_data.domain import summarize

class MaterialTests(unittest.TestCase):
    def test_primary_material(self):
        original,digest=load('data/normal.json');self.assertEqual(len(digest),64)
        value=summarize(original)
        self.assertEqual(value,summarize(json.loads(json.dumps(original))))
        self.assertEqual(digest,load('data/normal.json')[1])
        self.assertEqual(value['observed_total'],10);self.assertEqual(value['undirected_edges'],4)
    def test_invalid_material(self):
        for p in Path('data').glob('invalid-*.json'):
            with self.subTest(file=p.name):
                with self.assertRaises((ValueError,KeyError,TypeError)):summarize(load(p)[0])
    def test_additional_material(self):
        for p in Path('data').glob('*.json'):
            if p.name.startswith('invalid-'):continue
            with self.subTest(file=p.name):self.assertIsNotNone(summarize(load(p)[0]))
    def test_numeric_contract(self):
        for x in (True,float('inf'),float('nan'),'3'):
            with self.assertRaises(ValueError):finite(x)
        for x in (True,1.2,-1):
            with self.assertRaises(ValueError):integer(x)
        self.assertEqual(str(rational('2/6')),'1/3')
        for x in ('1/0',True,0.5):
            with self.assertRaises(ValueError):rational(x)
    def test_invalid_identity(self):
        for rows in ([{'id':'same'},{'id':'same'}],[{'id':''}]):
            with self.assertRaises(ValueError):unique(rows,'id')
    def test_json_rejection_and_utf8(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'资料.json';p.write_text('{"中心":"合成"}',encoding='utf8')
            self.assertEqual(load(p)[0],{'中心':'合成'})
            for raw in ('[]','{"a":NaN}','{"a":Infinity}'):
                p.write_text(raw,encoding='utf8')
                with self.assertRaises(ValueError):load(p)
