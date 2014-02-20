const char* dir = "/pool/PROOF-AAF/proof/dataset/MUON/ppillot";
//const char* dir = "/pool/PROOF-AAF/proof/dataset/default/cynthia";

void perm()
{
	Int_t n(1);

	while (n<10000)
	{  
		gProof->Exec(Form(".!chmod 777 %s",dir));
		gProof->Exec(Form(".!ls -ald %s",dir));
		gSystem->Sleep(10000);
		std::cout << "Loop " << n << std::endl;
		gProof->Exec(Form(".!ls -ald %s",dir));
		++n;
	}
}
	
