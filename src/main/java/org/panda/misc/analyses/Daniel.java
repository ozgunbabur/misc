package org.panda.misc.analyses;

import org.panda.resource.proteomics.RPPAIDMapper;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class Daniel
{
	public static final String BASE = "/Users/ozgun/Documents/Analyses/DanielDerrick/GPX4_casualPathInput_derrick/";

	public static void main(String[] args) throws IOException
	{
//		reformatData();
		prepareFolders();
	}

	private static void reformatData() throws IOException
	{
		RPPAIDMapper.convertFromMDACCToCP(BASE + "ACHN_GPX4_RPPA_limmaPvals_vsKOW_Derrick.tsv", BASE + "data2.txt");
	}

	private static void prepareFolders() throws IOException
	{
		String[] header = Files.lines(Paths.get(BASE + "data.txt")).findFirst().get().split("\t");

		for (int i = 4; i < header.length; i++)
		{
			String dirName = BASE + header[i];
			Files.createDirectory(Paths.get(dirName));

			BufferedWriter writer = Files.newBufferedWriter(Paths.get(dirName + "/parameters.txt"));

			Files.lines(Paths.get(BASE + "parameters.txt.template")).forEach(l -> FileUtil.writeln(l, writer));

			writer.write("value-column = " + header[i]);

			writer.close();
		}
	}
}
