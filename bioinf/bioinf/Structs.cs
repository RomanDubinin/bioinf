using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;

namespace bioinf
{
	public class TransitionDataColumn
	{
		public List<char> Match;
		public List<List<char>> Insertion;

		public TransitionDataColumn()
		{
			Match = new List<char>();
			Insertion = new List<List<char>>();
		}
	}

	public static class ExtensionMethods
	{
		// Deep clone
		public static T DeepClone<T>(this T a)
		{
			using (MemoryStream stream = new MemoryStream())
			{
				BinaryFormatter formatter = new BinaryFormatter();
				formatter.Serialize(stream, a);
				stream.Position = 0;
				return (T)formatter.Deserialize(stream);
			}
		}
	}
}