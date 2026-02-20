void skim() {
    // Open the original file and get the tree
    TFile *oldfile = new TFile("10x130ep_00mrad/lAger_v3.6.1_DVMP_JPsi_10x130ep_q2_1to50.hepmc3.tree.root");
    TTree *oldtree = (TTree*)oldfile->Get("hepmc3_tree");
    
    if (!oldtree) {
        std::cerr << "Error: Could not find 'hepmc3_tree' in the file." << std::endl;
        return;
    }

    Long64_t nentries = oldtree->GetEntries();

    // Create a new file and clone the tree structure (but keep 0 entries initially)
    TFile *newfile = new TFile("10x130ep_00mrad/lAger_v3.6.1_DVMP_JPsi_1ifb_10x130ep_q2_1to50.hepmc3.tree.root", "recreate");
    TTree *newtree = oldtree->CloneTree(0);

    // Loop over the events and only fill every 10th event
    for (Long64_t i = 0; i < nentries; i++) {
        if (i % 10 == 0) {
            oldtree->GetEntry(i);
            newtree->Fill();
        }
    }

    // Save and close
    newtree->AutoSave();
    delete oldfile;
    delete newfile;
    
    std::cout << "Skimming complete. Output saved to 10x130ep_00mrad/lAger_v3.6.1_DVMP_JPsi_1ifb_10x130ep_q2_1to50.hepmc3.tree.root" << std::endl;
}