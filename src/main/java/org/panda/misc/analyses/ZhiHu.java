package org.panda.misc.analyses;

import org.panda.resource.proteomics.MDACCFormatRPPALoader;
import org.panda.resource.proteomics.RPPAIDMapper;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

public class ZhiHu
{
	public static final String IN_BASE = "/home/ozgun/Data/ZhiHu-prostate-BRD4/";
	public static final String OUT_BASE = "/home/ozgun/Analyses/ZhiHu-prostate-BRD4/";

	public static void main(String[] args) throws IOException
	{
//		reformatData();
		prepareFolders();
	}

	private static void reformatData() throws IOException
	{
		RPPAIDMapper.convertFromMDACCToCP(IN_BASE + "PC3_JQ1_RPPA_limma_pvals.tsv", OUT_BASE + "data.txt");
	}

	private static void prepareFolders() throws IOException
	{
		String[] header = Files.lines(Paths.get(OUT_BASE + "data.txt")).findFirst().get().split("\t");

		for (int i = 4; i < header.length; i++)
		{
			String dirName = OUT_BASE + header[i];
			Files.createDirectory(Paths.get(dirName));

			BufferedWriter writer = Files.newBufferedWriter(Paths.get(dirName + "/parameters.txt"));

			Files.lines(Paths.get(OUT_BASE + "parameters.txt.template")).forEach(l -> FileUtil.writeln(l, writer));

			writer.write("value-column = " + header[i]);

			writer.close();
		}
	}
}
