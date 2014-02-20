/*
 *  ConfigureCuts.C
 *  aliroot
 *
 *  Created by Philippe Pillot on 01/04/10.
 *  Copyright 2010 SUBATECH. All rights reserved.
 *
 */

void ConfigureCuts(AliRunTagCuts *runCuts,
                   AliLHCTagCuts *lhcCuts,
                   AliDetectorTagCuts *detCuts,
                   AliEventTagCuts *evCuts)
{
  // Configure cuts.
  evCuts->SetNFWMuonRange(1,999999);
}

