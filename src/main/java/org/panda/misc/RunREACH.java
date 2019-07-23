package org.panda.misc;

import org.clulab.odin.Mention;
import org.clulab.processors.server.ProcessorServer;
import org.clulab.reach.PaperReader;
import org.clulab.reach.RunReachCLI;
import scala.collection.Iterator;
import scala.collection.Seq;
import scala.collection.immutable.List;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.zip.GZIPInputStream;

public class RunREACH
{
	public static void main(String[] args) throws IOException
	{
		deleteSmallFiles("/home/ozgun/Documents/reach/frext");
//		extractAbstracts("/media/ozgun/6TB/pubmed-abstracts-raw/pubmed19n0827.xml.gz", "/home/ozgun/Documents/reach/papers");
	}

	private static void deleteSmallFiles(String dir)
	{
		for (File file : new File(dir).listFiles())
		{
			if (file.length() < 100) file.delete();
		}
	}

	private static void extractAbstracts(String inFile, String outDir) throws IOException
	{
		BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inFile))));

		String line = reader.readLine();

		String pmid = null;

		while (line != null)
		{
			if (line.startsWith("      <PMID"))
			{
				pmid = getValue(line);
			}
			else if (line.equals("        <Abstract>"))
			{
				if (pmid == null)
				{
					System.err.println("An abstract has no PMID: " + reader.readLine());
				}
				else
				{
					line = reader.readLine();
					StringBuilder ab = new StringBuilder();
					while (!line.equals("        </Abstract>"))
					{
						if (line.startsWith("          <AbstractText>"))
						{
							ab.append(getValue(line));
						}
						else if (line.startsWith("          <AbstractText"))
						{
							String pre = line.substring(0, line.indexOf(">"));
							if (pre.contains("RESULTS") || pre.contains("CONCLUSIONS"))
							{
								ab.append(" ").append(getValue(line));
							}
						}
						line = reader.readLine();
					}
					String abs = ab.toString().trim();

					if (abs.isEmpty())
					{
						System.out.println(pmid + " has no abstract.");
					}
					else
					{
						BufferedWriter writer = Files.newBufferedWriter(Paths.get(outDir + "/" + pmid + ".txt"));
						writer.write(abs);
						writer.close();
					}
				}

				pmid = null;
			}

			line = reader.readLine();
		}

		reader.close();
	}

	private static String getValue(String line)
	{
		int start = line.indexOf(">");
		int end = line.indexOf("</");

		if (start < 0 || end < 0)
		{
			System.out.println("Malformed line = " + line);
			return "";
		}

		return line.substring(start + 1, end);
	}
}
