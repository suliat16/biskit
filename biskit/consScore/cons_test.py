#!/usr/bin/env python3
"""
"""
import os
from consScore import aminoCons as am
import biskit.test
from biskit import tools
from consScore import constool
from consScore import oma
from requests import exceptions
from unittest.mock import patch, MagicMock
from consScore import seq2conservation


class test_amino_conservation(biskit.test.BiskitTest):
    """
    Test suite testing the behaviour of the aminoCons module
    """

    TAGS = [biskit.test.EXE, biskit.test.LONG]

    @classmethod
    def setUpClass(cls):
        cls.filepath = tools.testRoot() + "consScore"

    def test_clean_alignment(self):
        """Tests to see that alignment file and the folder created are deleted
        after calling clean argument"""
        extra_aln = am.build_alignment(self.filepath + os.sep + 'multiFasta.fasta')
        self.assertTrue(os.path.exists(extra_aln))
        am.clean_alignment(extra_aln, cache=False)
        self.assertFalse(os.path.exists(extra_aln))

    def test_exe_mafft(self):
        """Tests that the Mafft executor module calls mafft"""
        extra_aln = am.MAFFT(self.filepath + os.sep + 'multiFasta.fasta')
        aln_path = extra_aln.run()
        self.assertTrue(os.path.exists(aln_path))
        am.clean_alignment(aln_path, cache=False)
        self.assertFalse(os.path.exists(aln_path))

    def test_get_alpha_g(self):
        """Tests that get alpha can retrieve the alpha parameter from the file"""
        digit = am.get_alpha(self.filepath + os.sep + 'multiFasta.res')
        self.assertEqual(digit, 2.83688)

    def test_get_alpha_nofile(self):
        """Test that an error is raised if a non-existent file is given to get_alpha"""
        with self.assertRaises(FileNotFoundError):
            am.get_alpha(self.filepath + os.sep + 'TheTree.txt')

    def test_getalpha_badfile(self):
        """Tests that an error is raised if get_alpha is called on an incompatible file"""
        with self.assertRaises(am.Rate4SiteError) as cm:
            am.get_alpha(self.filepath + os.sep + 'Fak2Human.fasta')
        err = cm.exception
        self.assertEqual(str(err), 'File format is not supported')

    def test_r2mat_g(self):
        """Tests that r2mat outputs the correct defaults for the dictionary"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta')
        output = test.rate2dict(self.filepath + os.sep + 'multiFasta.res')
        dictionary = {0: ('A', 0.6979),
                      1: ('A', 0.6979),
                      2: ('C', 0.7026),
                      3: ('C', 0.7026),
                      4: ('G', 0.1769),
                      5: ('G', 0.1769),
                      6: ('T', -1.577),
                      7: ('T', -1.577)}
        self.assertDictEqual(dictionary, output)

    def test_r2mat_multi(self):
        """Tests that the correct dictionary is output when some parameters are set"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta', score=True, qqint=True, gapped=True)
        output = test.rate2dict(self.filepath + os.sep + 'multiFasta.res')
        dictionary = {0: ('A', 0.6979, (-1.946, 2.735), '3/3'),
                      1: ('A', 0.6979, (-1.946, 2.735), '3/3'),
                      2: ('C', 0.7026, (-1.946, 2.735), '3/3'),
                      3: ('C', 0.7026, (-1.946, 2.735), '3/3'),
                      4: ('G', 0.1769, (-2.332, 1.853), '3/3'),
                      5: ('G', 0.1769, (-2.332, 1.853), '3/3'),
                      6: ('T', -1.577, (-3.889, -0.7852), '3/3'),
                      7: ('T', -1.577, (-3.889, -0.7852), '3/3')}
        self.assertDictEqual(output, dictionary)

    def test_r2mat_score(self):
        """Tests that correct dictionary is output when other parameters are set"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta', score=False, qqint=True, gapped=True, std=True)
        output = test.rate2dict(self.filepath + os.sep + 'multiFasta.res')
        dictionary = {0: ('A', (-1.946, 2.735), 2.836, '3/3'),
                      1: ('A', (-1.946, 2.735), 2.836, '3/3'),
                      2: ('C', (-1.946, 2.735), 2.836, '3/3'),
                      3: ('C', (-1.946, 2.735), 2.836, '3/3'),
                      4: ('G', (-2.332, 1.853), 2.725, '3/3'),
                      5: ('G', (-2.332, 1.853), 2.725, '3/3'),
                      6: ('T', (-3.889, -0.7852), 2.309, '3/3'),
                      7: ('T', (-3.889, -0.7852), 2.309, '3/3')}
        self.assertDictEqual(output, dictionary)

    def test_r2prof(self):
        """tests that rate2profile correctly creates profile collections with the correct output"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta', qqint=True, std=True, gapped=True)
        output = test.rate2profile(self.filepath + os.sep + 'multiFasta.res')
        self.assertEqual(5, len(output))
        self.assertTrue('Amino Acid' in output)
        self.assertTrue('Conservation Score' in output)
        self.assertTrue('Standard Deviation' in output)
        self.assertTrue('QQ interval' in output)
        self.assertTrue('Gapped' in output)
        self.assertTrue(isinstance(output, biskit.ProfileCollection))

    def test_r2prof_default(self):
        """Tests that the default rate2site instantiation produces the right profile"""
        test = am.Rate4Site(self.filepath + os.sep + 'multiFasta.fasta')
        output = test.rate2profile(self.filepath + os.sep + 'multiFasta.res')
        self.assertTrue('Amino Acid' in output)
        self.assertTrue('Conservation Score' in output)
        self.assertEqual(2, len(output))
        self.assertTrue(isinstance(output, biskit.ProfileCollection))

    def test_r4s_close(self):
        """Tests to see that close deletes the correct files"""
        r4sobject = am.Rate4Site(self.filepath + os.sep + 'multiFasta.aln')
        r4sobject.run()
        self.assertTrue(os.path.isfile(os.getcwd() + os.sep + 'multiFasta.res'))
        r4sobject.close()
        self.assertFalse(os.path.isfile(os.getcwd() + os.sep + 'multiFasta.res'))
        self.assertFalse(os.path.isdir(os.getcwd() + os.sep + 'multiFasta'))

    @classmethod
    def tearDownClass(cls):
        am.clean_alignment(os.getcwd() + os.sep + 'multiFasta.aln', cache=False)


class test_constool(biskit.test.BiskitTest):

    """Test suite testing the behaviour of the constool module"""

    TAGS = [biskit.test.NORMAL]

    def test_get_num_good(self):
        """Tests that get num can retrieve a single float from a string"""
        digit = constool.get_num("Hello, how are all 1234.56 of you today")
        self.assertEqual(digit[0], 1234.56)

    def test_get_num_multi(self):
        """Tests that get_num can retrieve multiple floats and integers from a string """
        digit = constool.get_num("Hello, how are all 1234.56 and 42 of you today")
        self.assertEqual(digit, [1234.56, 42])

    def test_get_num_none(self):
        """Tests that a string containing no digits returns no digits"""
        digit = constool.get_num("Oh, there are none of you today")
        self.assertEqual(digit, [])

    def test_get_num_negative(self):
        """Tests that get_num can retrieve negative and positive numbers from within a string"""
        digit = constool.get_num("Hello, how are all 1234.56 and -42 of you today")
        self.assertEqual(digit, [1234.56, -42])

    def test_header_chk(self):
        """test that header_check adds a header to a header-free string """
        output = constool.header_check("Hello World")
        self.assertEqual(output, ">Input Sequence\nHello World")

    def test_ortho_empty(self):
        """tests that header_check raises an exception if an empty sequence is entered"""
        with self.assertRaises(constool.SequenceError) as cm:
            constool.header_check("")
        err = cm.exception
        self.assertEqual(str(err), "Empty Sequence entered.")

    def test_non_sequence(self):
        """tests that header_check raises an exception if given a non alphabetic sequence"""
        with self.assertRaises(constool.SequenceError) as sm:
            constool.header_check("003893")
        err = sm.exception
        self.assertEqual(str(err), "Not a sequence. Please try again")

    def test_ortho_already_fasta(self):
        """tests that header_check doesn't change fasta sequences that already have a header"""
        output = constool.header_check(">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")
        self.assertEqual(output, ">Pingo Pongo | Happiness\nWOEFJEKTJEJTEK")

    def test_get_fasta_seq(self):
        """Tests that get_fasta_seq only gets the sequence of the fasta string, not the identifying line"""
        tester = constool.get_fasta_sequence(""">OAP01791.1 CDC48A [Arabidopsis thaliana]
                                    MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQLFRGDTILIKGKKRKD
                                    TVCIALADETCEEPKIRMNKVVRSNLRVRLGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLK""")
        self.assertFalse('>OAP01791.1 CDC48A [Arabidopsis thaliana]' in tester)
        self.assertTrue('SKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQL' in tester)

    def test_get_fasta_multi(self):
        """tests that get_fasta returns a list of sequences given a fasta file"""
        tester = constool.get_fasta_sequence(""">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]
MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE
>LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]
MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEASTPKVSKQGRSEEISESE
>ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]
MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQSAKKARVEEASTPKVNKQSRSEXETSAP""", index=1)
        self.assertFalse("[Loxodonta africana]" in tester)
        self.assertFalse(">PROCA12070 | ENSPCAG00000012030" in tester)
        self.assertEqual(tester, "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEASTPKVSKQGRSEEISESE")

    def test_indv_blk(self):
        """tests that indv_block can pull out individual fasta sequences with their identifying line"""
        tester = constool.indv_block(""">OAP01791.1 CDC48A [Arabidopsis thaliana]
                                    MSTPAESSDSKSKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQLFRGDTILIKGKKRKD
                                    TVCIALADETCEEPKIRMNKVVRSNLRVRLGDVISVHQCPDVKYGKRVHILPVDDTVEGVTGNLFDAYLK""")
        self.assertTrue('>OAP01791.1 CDC48A [Arabidopsis thaliana]' in tester[0])
        self.assertTrue('SKKDFSTAILERKKSPNRLVVDEAINDDNSVVSLHPATMEKLQL' in tester[0])

    def test_indv_blk_multi(self):
        """tests that indv_block can pull out individual fasta sequences from a string with multiple"""
        tester = constool.indv_block(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
"MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
"MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"  
"MKTRQNKDSMSMRSGRKKEAPGPREELRS")
        self.assertEqual(len(tester), 3)
        self.assertEqual(tester[1],(">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                    "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA" ))
    def test_remove_protein(self):
        """tests that remove_protein can accurately remove a  protein with the given id"""
        tester = constool.remove_protein(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS", "LOXAF14113")
        self.assertEqual(tester, ">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS")

    def test_remove_protein_nothing(self):
        """Tests that remove protein deosnt remove anythign if the entry ID is not in string"""
        tester = constool.remove_protein(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS", "PROXY90210")
        self.assertEqual(tester, ">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS")

    def test_remove_first_protein(self):
        """Tests that remove_first_protein removes the first protein in a fasta format file"""
        tester = constool.remove_first_protein(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n"
                                     ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS")
        self.assertEqual(tester, ">LOXAF14113 | G3TAL7 | HOG:0377891.2a.2a | [Loxodonta africana]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTA\n"
                                     ">ECHTE02547 | ENSETEG00000016682 | HOG:0377891.2a.2a | [Echinops telfairi]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRS")

    def test_seqnwl_strip(self):
        """Tests that seqnwl_strip removes the newlines from within the fasta string, but keeps the newline at the end"""
        tester = constool.seqnwl_strip(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE\n")
        self.assertEqual(tester, ">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE")

    def test_seqnwl_strip_messy(self):
        """Tests that seqnwl_strip removes the newlines from within the fasta string"""
        tester = constool.seqnwl_strip(">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNK\nDSMSMRSGRKKEAPGPREEL\nRSRGRASPGGVSTSSSDGKAEKSRQTAK\nKARVEEVSAPKVSKQGRGEEIS\nESE\n")
        self.assertEqual(tester, ">PROCA12070 | ENSPCAG00000012030 | HOG:0377891.2a.2a | [Procavia capensis]\n"
                                     "MKTRQNKDSMSMRSGRKKEAPGPREELRSRGRASPGGVSTSSSDGKAEKSRQTAKKARVEEVSAPKVSKQGRGEEISESE")


class TestCons(biskit.test.BiskitTest):
    """
    Test suite testing the behaviour of the oma module
    """

    TAGS = [biskit.test.NORMAL]

    def setUp(self):
        lysozyme = ("MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGY"
                    "NTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVAC"
                    "AKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV")
        self.lyz = oma.OrthologFinder(lysozyme)
        stragg = ('YYYYYYYYYYYVVVVVVVVVVVVBBBBBBBBBBBAAAAAAAAAAANNNNNNNNNNN'
                  'MMMMMMMMMMMWWWWWWWWWWWWWCCCCCCCCCCSSSSSSSSSS')
        self.aggregate = oma.OrthologFinder(stragg)

    @classmethod
    def setUpClass(cls):
        cls.filepath = tools.testRoot() + "consScore"

        with open((cls.filepath + os.sep + 'CDC48Aseq.txt'), 'r') as file:
            arabidopsisCDC48A = file.read()
        cls.CDC48A = oma.OrthologFinder(arabidopsisCDC48A)

        with open((cls.filepath + os.sep + 'OMA_get_id.txt'), 'r') as file:
            OMA_id_response = file.read()
        cls.response = bytes(OMA_id_response, 'utf-8')

        with open((cls.filepath + os.sep + 'orhtolog_response.txt'), 'r') as file:
            ortho_response = file.read()
        cls.oresponse = bytes(ortho_response, 'utf-8')

        with open((cls.filepath + os.sep + 'fasta_response.txt'), 'r') as file:
            fasta_response = file.read()
        cls.fresponse = bytes(fasta_response, 'utf-8')

        with open((cls.filepath + os.sep + 'ATN1_HOGS.txt'), 'r') as file:
            HOG_response = file.read()
        cls.hresponse = bytes(HOG_response, 'utf-8')

        with open((cls.filepath + os.sep + 'hog_content.txt'), 'r') as file:
            level_response = file.read()
        cls.lvlresponse = bytes(level_response, 'utf-8')

    def test_read_resp_retOMA(self):
        """tests that read_resp_protID correctly parses the JSON response and retrieves the protein id"""
        the_response = MagicMock(content=self.response)
        self.CDC48A.read_resp_protID(the_response)
        self.assertTrue(self.CDC48A.id == "ARATH09528")

    def test_read_resp_upOMA(self):
        """tests that read_resp_protID correctly parses the JSON response and saves a list of the ortholog ids"""
        the_response = MagicMock(content=self.oresponse)
        self.CDC48A.read_resp_orthoIDs(the_response)
        self.assertTrue("H3DFZ9" in self.CDC48A.ortholog_ids)

    def test_build_url_retOMA(self):
        """Tests that build URL inserts the variable portions of the url correcly into the base"""
        self.lyz.sequence = 'MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV'
        self.assertEqual(constool.build_url(base_url='https://omabrowser.org', tail='/api/sequence/?query={0}',
                                            variation=[self.lyz.sequence]),
                         'https://omabrowser.org/api/sequence/?query=MKALIVLGLVLLSVTVQGKVFERCELARTLKRLGMDGYRGISLANWMCLAKWESGYNTRATNYNAGDRSTDYGIFQINSRYWCNDGKTPGAVNACHLSCSALLQDNIADAVACAKRVVRDPQGIRAWVAWRNRCQNRDVRQYVQGCGV')

    @patch('oma.requests.get')
    def test_ologs_stragg(self, requests_mock):
        """Tests that a request with a bad status code raises an exception with call_orthologs"""
        requests_mock.requests.get.return_value = None
        requests_mock.requests.get.status_code = 400
        with self.assertRaises(exceptions.RequestException) as cm:
            straggFasta = self.aggregate.get_orthologs()
        err = cm.exception
        self.assertTrue('There was an issue querying the database. Status code' in str(err))

    @patch('oma.requests.get')
    def test_ofasta_stragg(self, requests_mock):
        """Tests that a bad request raises an exception in ortholog_to_fasta"""
        requests_mock.requests.get.status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.ID_to_fasta()

    @patch('oma.requests.get')
    def test_ofasta_cdc(self, requests_mock):
        """Tests that ortholog_to_fasta correctly parses the response data"""
        requests_mock().status_code = 200
        requests_mock().text = self.fresponse
        test = self.CDC48A.ID_to_fasta()
        self.assertTrue('[Arabis alpina]' in test)

    @patch('oma.requests.get')
    def test_orIDs_stragg(self, requests_mock):
        """Test that update_orthoIDs returns an exception given a bad request status"""
        requests_mock().status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.update_orthoIDs()

    @patch('oma.requests.get')
    def test_ret_hogs_stragg(self, requests_mock):
        """Tests that retrieve_HOG_level throws an exception when given a bad request status"""
        requests_mock().status_code = 400
        with self.assertRaises(exceptions.RequestException):
            self.aggregate.retrieve_HOG_level()

    def test_empty_input(self):
        """Checks that the correct exception is raised when an empty sequence is entered"""
        fs = oma.OrthologFinder("")
        with self.assertRaises(oma.SequenceError) as cm:
            fs.get_orthologs()
        err = cm.exception
        self.assertEqual(str(err), 'Input sequence is empty!')

    def test_empty_input_HOG(self):
        """Checks that the correct exception is raised when an empty sequence is entered"""
        fs = oma.OrthologFinder("")
        with self.assertRaises(oma.SequenceError) as cm:
            fs.get_HOGs()
        err = cm.exception
        self.assertEqual(str(err), 'Input sequence is empty!')

    @patch('oma.requests.get')
    def test_hog_fasta(self, mock_request):
        """Tests that HOG_to_fasta correctly parses the request response"""
        self.lyz.hog_level = "Amniota"
        mock_request().status_code = 200
        mock_request().text = self.hresponse
        test = self.lyz.HOG_to_fasta()
        self.assertTrue("HOG:0377891.2a.2a" in test)

    @patch('oma.requests.get')
    def test_read_hog_roottrue(self, mock_request):
        """tests that read_HOGid retrieves the root ID when root=True"""
        thing = MagicMock(content=self.lvlresponse)
        mock_request().content = self.lvlresponse
        test = self.lyz.read_HOGid(thing, root=True)
        self.assertEqual(test, 'Amniota')

    @patch('oma.requests.get')
    def test_read_hog_rootfalse(self, mock_request):
        """tests that read_HOGid retrieves a list of the ids when root=False """
        thing = MagicMock(content=self.lvlresponse)
        mock_request().content = self.lvlresponse
        test = self.lyz.read_HOGid(thing, root=False)
        self.assertTrue('Hystricomorpha' in test)
        self.assertTrue('Caniformia' in test)
        self.assertTrue('Gorilla gorilla gorilla' in test)

class Test(biskit.test.BiskitTest):
    """
    Test suite for the seq2conservation pipeline
    """

    TAGS = [biskit.test.NORMAL, biskit.test.LONG]

    @classmethod
    def setUpClass(cls):
        cls.filepath = tools.testRoot() + "consScore"

        with open((cls.filepath+os.sep+'CDC48Aseq.txt'), 'r') as file:
            arabidopsisCDC48A = file.read()
        cls.CDC48A = seq2conservation.ConservationPipe(arabidopsisCDC48A, cache=False)

        with open((cls.filepath+os.sep+'multiFasta.fasta'), 'r') as file:
            cls.ex_seq = file.read()

    @patch('seq2conservation.oma.OrthologFinder.get_HOGs')
    def test_call_orthologs(self, HOG_mock):
        """tests that call_orthologs correctly calls methods to make an output file of 'orthologs'"""
        HOG_mock.return_value = self.ex_seq
        tester = self.CDC48A.call_orthologs()
        self.assertTrue(os.path.isfile(tester))
        self.assertTrue('Protein_Sequence.orth' in tester)

    @patch('seq2conservation.aminoCons.build_alignment')
    def test_call_alignment(self, mock_aln):
        """tests that call_alignment calls the correct methods and generates the correct output"""
        mock_aln.return_value = self.filepath+os.sep + 'multiFasta.aln'
        tester = self.CDC48A.call_alignment(self.filepath+os.sep+'multiFasta.fasta')
        self.assertTrue(mock_aln.called)
        self.assertTrue(os.path.isfile(tester))
        self.assertTrue('multiFasta.aln' in tester)
        am.clean_alignment(tester, cache=False)

    @patch('seq2conservation.aminoCons.get_alpha')
    def test_call_rate4site(self, mock_score):
        """tests that call_rate4site calls the correct methods and generates the correct output"""
        mock_score.return_value = 2.83688
        tester= self.CDC48A.call_rate4site(os.getcwd()+os.sep+'example_data'+os.sep + 'multiFasta.aln')
        self.assertTrue(mock_score.called)
        self.assertEqual(2.83688, tester)

    @patch('seq2conservation.oma.OrthologFinder.get_HOGs')
    @patch('seq2conservation.aminoCons.build_alignment')
    @patch('seq2conservation.aminoCons.Rate4Site.run')
    def test_pipe(self, mock_run, mock_aln, mock_hog):
        """tests that the pipe calls the correct methods and generates the correct output"""
        mock_hog.return_value = self.ex_seq
        mock_aln.return_value = os.getcwd()+os.sep+'example_data'+os.sep + 'multiFasta.aln'
        mock_run.return_value = {0:('A', 0.6979),
                    1:('A', 0.6979),
                    2:('C', 0.7026),
                    3:('C', 0.7026),
                    4:('G', 0.1769),
                    5:('G', 0.1769),
                    6:('T', -1.577),
                    7:('T', -1.577)}
        tester = self.CDC48A.pipe()
        self.assertTrue(mock_run.called)
        self.assertTrue(mock_aln.called)
        self.assertTrue(mock_hog.called)
        self.assertTrue(type(tester), dict)

    @classmethod
    def tearDownClass(cls):
        os.remove(os.getcwd()+ os.sep + 'Protein_Sequence.orth')


if __name__ == '__main__':
    biskit.test.localTest()
