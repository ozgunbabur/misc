package org.panda.misc;

import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * For various conversions.
 * @author Ozgun Babur
 */
public class CausalPathFormatConvertor
{
	public static void main(String[] args) throws IOException
	{
		convertPHOTONToCP("/home/babur/Documents/Analyses/CausalPathAlternatives/PHOTON/EGFData/data.csv",
			"/home/babur/Documents/Analyses/CausalPathAlternatives/PHOTON/EGFData/CausalPath-analysis/data.txt");
	}
	public static void convertPHOTONToCP(String inFile, String outFile) throws IOException
	{
		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));
		writer.write("ID\tSymbols\tSites\tEffect\tValue");
		Files.lines(Paths.get(inFile)).skip(1).map(l -> l.split(",")).forEach(t ->
			FileUtil.lnwrite(t[4] + "-" + t[1] + t[2] + "\t" + t[4] + "\t" + t[1] + t[2] + "\t\t" + t[3], writer));
		writer.close();
	}
}
