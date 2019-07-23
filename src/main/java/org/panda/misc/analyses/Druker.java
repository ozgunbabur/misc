package org.panda.misc.analyses;

import org.panda.resource.proteomics.MDACCFormatRPPALoader;
import org.panda.resource.proteomics.RPPAIDMapper;

import java.io.IOException;

public class Druker
{
	public static void main(String[] args) throws IOException
	{
		String dir = "/home/ozgun/Analyses/Druker/Convert-from-MDACC/";
		convertFromMDACC(dir + "20_Jeff_Tyner__Anna_Schultz_C_normlogmcentered.txt", dir + "dataC.txt");
	}

	public static void convertFromMDACC(String inFile, String outFile) throws IOException
	{
		RPPAIDMapper.convertFromMDACCToCP(inFile, outFile);
	}
}
