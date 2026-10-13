"""Negative checks: corrupted inputs must not produce a valid certificate.

Run: sage -python F2_human_audit/test_certificate.py
"""
import json
from copy import copy
from pathlib import Path
import tempfile
import unittest

import sage.all as sage
import certify


class CertificateRejectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.audit = certify.load_namespace('check_F2_product_properties.sage')

    def test_changed_printed_matrix_is_rejected(self):
        changed = dict(self.audit)
        changed['extension_matrix_R'] = copy(self.audit['extension_matrix_R'])
        changed['extension_matrix_R'][64, 0] += 1
        with self.assertRaisesRegex(RuntimeError, 'printed paper matrix'):
            certify.check_paper_data(changed)

    def test_changed_paper_eta_is_rejected(self):
        snapshot = json.loads((certify.AUDIT / 'paper_instance.json').read_text())
        text = '\n'.join(f'{i} & \\texttt{{{v}}}' for i, v in enumerate(snapshot['R_rows_u64']))
        text += '\n\\label{eq:eta}\n\\eta(x)=x^2+x.\n\\end{equation}'
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'explicit_instance.tex'
            path.write_text(text)
            with self.assertRaisesRegex(RuntimeError, 'supplied paper eta'):
                certify.check_paper_data(self.audit, path)

    def test_nonsystematic_matrix_is_rejected(self):
        changed = dict(self.audit)
        changed['extension_matrix_R'] = copy(self.audit['extension_matrix_R'])
        changed['extension_matrix_R'][0, 1] = 1
        with self.assertRaisesRegex(RuntimeError, 'not systematic'):
            certify.check_structure(changed)

    def test_changed_product_prefix_is_rejected(self):
        changed = dict(self.audit)
        original = self.audit['product_rule']
        def corrupt(left, right):
            result = original(left, right)
            result[0] += left[1]*right[1]
            return result
        changed['product_rule'] = corrupt
        with self.assertRaisesRegex(RuntimeError, 'block convolution mismatch'):
            certify.check_structure(changed)

    def test_duplicate_evaluator_is_rejected(self):
        # The function closes over the module namespace; restore it even if a
        # check fails so another test never inherits the corruption.
        original = self.audit['gamma_numerators']
        changed = list(original)
        changed[1] = changed[0]
        self.audit['gamma_numerators'] = changed
        try:
            with self.assertRaisesRegex(RuntimeError, 'numerator rank'):
                self.audit['check_nondegeneracy'](show_progress=False)
        finally:
            self.audit['gamma_numerators'] = original

    def test_changed_direct_base_evaluator_is_rejected(self):
        changed = dict(self.audit)
        changed['base_gamma_numerators'] = [list(row) for row in self.audit['base_gamma_numerators']]
        changed['base_gamma_numerators'][0][0] += 1
        with self.assertRaisesRegex(RuntimeError, 'direct base evaluator differs'):
            certify.check_direct_base_evaluator(changed)

    def test_extension_map_preserves_binary_list_interface(self):
        result = self.audit['extend_base_message']([1]*64)
        self.assertEqual(result.base_ring(), sage.GF(2))
        self.assertEqual(list(result[:64]), [1]*64)


if __name__ == '__main__':
    unittest.main()
