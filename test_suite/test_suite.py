#!/usr/bin/env python
"""
 Test suite for the Quest program.

 The suite consists of a series of different types of tests. The simplest
 (and the majority of the tests) are energy tests where a quest input
 file is run via quest_run(), and the energy is compared to a reference
 value.

 There are also simple tests to see whether an MD run completed correctly.

 These runs are all submitted through the Python Unit Test framework.

 Written by Rick Muller
"""

import unittest,time,os
from quest import quest_run,quest_md_isok

class QuestSuite(unittest.TestCase):
    def testH2(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('h2.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-1.786568,4)
      
    def testH2Geopt(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('h2_geopt.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-1.781264,4)
        
    def testH2GGA(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('h2_gga.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-1.8399273,4)
        
    def testH2GeoptGGA(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('h2_geopt_gga.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-1.835298,4)
        
    def testSi(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('si_111.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-62.258784,4)
      
    def testSiIdosolv3(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('si_111_idosolv3.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-62.258784,4)
      
    def testSi112(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('si_112.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-125.257398,4)
      
    def testSiComplex(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('si_111_complex.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-62.269826,4)
      
    def testSiCell(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('si_111_cell.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-62.308584,4)
      
    def testGaN(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('gan.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-58.798295,4)

class LongQuestSuite(unittest.TestCase):
    def testBenzene(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('c6h6.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-75.420703,4)
      
    def testRu(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('ru_slab.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-232.249206,4)

    # Commented out, because the test is unreliable and takes a long time
    #def testNEB(self,**kwargs):
    #    kwargs['do_preclean'] = True
    #    kwargs['do_clean'] = True
    #    kwargs['do_setup'] = True
    #    energy = quest_run('hclh_neb.in',**kwargs)
    #    ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
    #    ofile.flush()
    #    self.assertAlmostEqual(energy,-31.996766,4)

    def testAs(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('as_test.in',**kwargs)
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-28.3216914,4)
        
    # Note that the tolerances on the MD runs are weaker
    def testMD(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        energy = quest_run('h2o-dimer-md.in',**kwargs)
        md_isok = quest_md_isok('h2o-dimer-md.in')
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-68.782839,3)
        self.assert_(md_isok)

    # Note that the tolerances on the MD runs are weaker
    def testMDRestart(self,**kwargs):
        kwargs['do_preclean'] = True
        kwargs['do_clean'] = True
        kwargs['do_setup'] = True
        kwargs['remove_vxyz'] = False
        energy = quest_run('h2o-dimer-md-restart.in',**kwargs)
        md_isok = quest_md_isok('h2o-dimer-md-restart.in')
        ofile.write("%-40s  %20.10f\n" % (self.id(),energy))
        ofile.flush()
        self.assertAlmostEqual(energy,-68.784460,3)
        self.assert_(md_isok)


if __name__ == '__main__':
    import getopt,sys
    global ofile

    opts,args = getopt.getopt(sys.argv[1:],'alsv')

    runlevel = 0
    verbosity = 1
    for key,val in opts:
        if key == '-v':
            verbosity = 2
        if key == '-s':
            runlevel = 0
        if key == '-l':
            runlevel = 1
        if key == '-a':
            runlevel = 1
    
    ofile = open("test-results.txt","w")
    host = os.uname()[1]
    ofile.write("Test suite run %s on %s\n\n" % (time.asctime(),host))
    ofile.flush()
    # Verbose version:
    shorttests = unittest.makeSuite(QuestSuite)
    longtests = unittest.makeSuite(LongQuestSuite)
    if runlevel == 0:
        suite = shorttests
    else:
        suite = unittest.TestSuite((shorttests,longtests))
    unittest.TextTestRunner(verbosity=verbosity).run(suite)
    # Non-verbose version:
    #unittest.main()
    ofile.close()
    
