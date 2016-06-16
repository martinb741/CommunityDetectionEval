package cdTest;

import java.util.concurrent.TimeUnit;

import javax.ws.rs.client.Client;
import javax.ws.rs.client.Entity;
import javax.ws.rs.client.WebTarget;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;

import org.glassfish.jersey.apache.connector.ApacheConnectorProvider;
import org.glassfish.jersey.client.ClientConfig;
import org.glassfish.jersey.client.ClientProperties;
import org.glassfish.jersey.client.JerseyClientBuilder;

import com.google.common.base.Stopwatch;

public class OCDREST {

	public static void main(String[] args) {
		
		int[] graphs = new int[] {107};
		String[] algorithms = new String[] {"SPEAKER_LISTENER_LABEL_PROPAGATION_ALGORITHM"};
		
		for (int graph : graphs){
			for (String algorithm : algorithms){
				String parameters = getParameters(algorithm);
				runOCD(graph, algorithm, parameters);
			}
		}
		
	}
	
	private static void runOCD(int graph, String algorithm, String parameters){
		/*
		 * Run an OCD algorithm on WebOCD
		 * 
		 * example curl command to run SLPA on graph number 7
		 * curl -X POST --header 'Content-Type: text/plain' --header 'Accept: text/xml'
		 *      -d '<Parameters><Parameter><Name>memorySize</Name><Value>100</Value></Parameter><Parameter><Name>probabilityThreshold</Name><Value>0.15</Value></Parameter></Parameters>'
		 *      'http://localhost:8080/ocd/covers/graphs/7/algorithms?name=SLPA_Cover&algorithm=SPEAKER_LISTENER_LABEL_PROPAGATION_ALGORITHM&contentWeighting=false&componentNodeCountFilter=0'
		 *
		 * Response:
		 * <Cover><Id><CoverId>3</CoverId><GraphId>7</GraphId></Id></Cover> 
		 */
		
        ClientConfig config = new ClientConfig();
        config.connectorProvider(new ApacheConnectorProvider());
        config.property(ClientProperties.PROXY_URI, "http://localhost:8082");
		Client client = JerseyClientBuilder.newClient(config);
		
		WebTarget target = client.target("http://localhost:8080")
				.path("ocd").path("covers").path("graphs").path(Integer.toString(graph)).path("algorithms")
				.queryParam("name", "TestCover")
				.queryParam("algorithm", algorithm)
				.queryParam("contentWeighting", "false")
				.queryParam("componentNodeCountFilter", "0");
		
		Stopwatch sw = Stopwatch.createStarted();

		Response response = target.request(MediaType.TEXT_XML).post(Entity.entity(parameters, MediaType.TEXT_PLAIN));
		
		System.out.println("Response Status: " + response.getStatus());
		System.out.println("Response Text: " + response.readEntity(String.class));
		
		sw.stop();
		
		long computeTime = sw.elapsed(TimeUnit.SECONDS);
		
		System.out.println("OCD algorithm: " + algorithm + " / computation time: " + computeTime);
	}
	
	private static String getParameters(String algorithm){
		
		String parameters = "";
		
		switch (algorithm){
		case "SPEAKER_LISTENER_LABEL_PROPAGATION_ALGORITHM":
			parameters = "<Parameters><Parameter><Name>memorySize</Name><Value>100</Value></Parameter><Parameter><Name>probabilityThreshold</Name><Value>0.15</Value></Parameter></Parameters>";
			break;
		}
		
		return parameters;
	}

}
