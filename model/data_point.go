package model

type DataPoint struct {
	Name            string   // Name of the place
	Lat             float64  // Latitude of the place
	Lng             float64  // Longitude of the place
	Categories      []string // A list of category the place is classified under
	AttributeVector []bool   // A boolean vector indicating if an attribute is true or false for this place. The size
	// of the attribute vector should be the same for all data model, since each index maps to a unique attribute.
}
