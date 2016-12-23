package org.panda.misc;

import org.biopax.paxtools.controller.PathAccessor;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.Complex;
import org.biopax.paxtools.model.level3.Control;
import org.biopax.paxtools.model.level3.Interaction;
import org.biopax.paxtools.model.level3.Process;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * For mining out new information from BioPAX, given a SIF file.
 * @author Ozgun Babur
 */
public class EnrichSIFFile
{
	private static final PathAccessor xrefAcc = new PathAccessor("XReferrable/xref:PublicationXref/id");
	private static final PathAccessor evidAcc = new PathAccessor("Observable/evidence/xref:PublicationXref/id");

	public static void main(String[] args) throws IOException
	{
		String base = "/home/babur/Documents/Analyses/MrOS/";
		addPubmedIDs(base + "paths-between-all.sif", base + "paths-between-all-pubmed.sif",
			"/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl");
	}

	public static void addPubmedIDs(String sifFile, String newFile, String biopaxFile) throws IOException
	{
		SimpleIOHandler io = new SimpleIOHandler();
		Model model = io.convertFromOWL(new FileInputStream(biopaxFile));

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(newFile));

		Files.lines(Paths.get(sifFile)).map(l -> l.split("\t")).forEach(t ->
		{
			FileUtil.write(t[0] + "\t" + t[1] + "\t" + t[2], writer);
			if (t.length > 3 && t[3] != null && !t[3].isEmpty())
			{
				List<String> list = new ArrayList<>(collectPubmedIDs(t[3].split(" |;"), model));
				Collections.sort(list);
				FileUtil.write("\t" + CollectionUtil.merge(list, ";"), writer);
			}
			FileUtil.write("\n", writer);
		});

		writer.close();
	}



	private static Set<String> collectPubmedIDs(String[] ids, Model model)
	{
		Set<String> set = new HashSet<>();
		for (String id : ids)
		{
			set.addAll(collectPubmedIDs(id, model));
		}
		return set;
	}

	private static Set<String> collectPubmedIDs(String id, Model model)
	{
		BioPAXElement ele = model.getByID(id);
		Set<String> set = new HashSet<>();

		if (ele instanceof Interaction)
		{
			set.addAll(xrefAcc.getValueFromBean(ele));
			set.addAll(evidAcc.getValueFromBean(ele));
		}

		if (ele instanceof Control)
		{
			for (Process pr : ((Control) ele).getControlled())
			{
				set.addAll(collectPubmedIDs(pr.getUri(), model));
			}
		}

		if (ele instanceof Complex)
		{
			for (Interaction inter : ((Complex) ele).getParticipantOf())
			{
				set.addAll(collectPubmedIDs(inter.getUri(), model));
			}
		}

		return set;
	}
}
